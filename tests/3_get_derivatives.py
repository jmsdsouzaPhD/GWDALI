import sys, os
import GWDALI as gw
import numpy as np
import pandas as pd
import jax.numpy as jnp
import matplotlib.pyplot as plt

from time import time as now
from tqdm import trange

np.random.seed(0)

Approxs = ["TaylorF2"] + [f"IMRPhenom{x}" for x in "A,B,C,D,HM".split(',')]

rad = np.pi/180
deg = 1./rad

detector = {'name':'ET','lon':0, 'lat':90, 'rot':0, 'shape':90}
detectors = [detector]

freq = 10**np.linspace(np.log10(5),3,1000)

Mc = 44.#36.
eta = 0.24
dL = 460.
iota = np.pi/3
m1 = 0.5*(Mc/eta**(3./5)) * (1 + np.sqrt(1-4*eta))
m2 = 0.5*(Mc/eta**(3./5)) * (1 - np.sqrt(1-4*eta))
q = m2/m1
M = m1+m2
sz1 = 0.5
sz2 = -0.2

AllMass = {"Mc": Mc, "eta": eta, "m1": m1, "m2": m2, "M": M, "q": q}
mass_keys = "Mc,eta,m1,m2,M,q".split(',')
n = len(mass_keys) ; k = 0
FreeParams_vec = []
for i in range(n):
	for j in range(i+1,n):
		fprms =  [ mass_keys[i], mass_keys[j] ]
		FreeParams_vec.append(fprms)
		#print(f"{k+1}\tFreeParams: {fprms}")
		k += 1

GwPrms = {}
GwPrms['dL']   = dL/1.e3 # Gpc
GwPrms['iota'] = iota 
GwPrms['psi']  = np.pi/5 
GwPrms['t_coal']   = 0. 
GwPrms['phi_coal'] = 0. 
GwPrms['RA']   = 45.
GwPrms['Dec']  = 45.
GwPrms["sz1"]  = sz1 
GwPrms["sz2"]  = sz2 

GwPrms["sx1"] = 0.
GwPrms["sy1"] = 0.
GwPrms["sx2"] = 0.
GwPrms["sy2"] = 0.

approx = "IMRPhenomHM"
diff_order  = "first"

#FreeParams = "Mc,eta".split(',')
#FreeParams = FreeParams_vec[int(sys.argv[1])-1]
#FreeParams = sys.argv[1].split(',')

#for m in FreeParams: GwPrms[m] = AllMass[m]

GwPrms["Mc"] = Mc ; GwPrms["eta"] = eta
FreeParams = "dL,iota,psi,phi_coal,RA,Dec,t_coal,Mc,eta,sz1,sz2".split(',')


#FreeParams = "dL,iota,psi".split(',')
#FreeParams = "dL,iota,psi,phi_coal,RA,Dec,t_coal,Mc,eta,sz1,sz2".split(',')
#FreeParams = "inv_dL,cos_iota,psi,phi_coal,RA,Dec,t_coal,Mc,inv_eta,sz1,sz2".split(',')

dets = [detector, detector]

fig = plt.figure(figsize=(8,10))
plt.rcParams.update({"font.family":"serif","font.size":20,"font.weight":"bold"})

#wftype = "lal"
'''
for wftype in ["lal","jax"]:
	if wftype == "lal":
		DiffMethods = ["numdiff"]
	else:
		DiffMethods = ["numdiff","autodiff"]

	for diff_method in DiffMethods:
'''

methods = [ ["jax","autodiff"], ["jax","numdiff"] , ["lal","numdiff"] ]
lss = "-","--",":"
clrs = "k","b","r"

methods.reverse()
if True:
	for ni, m in enumerate(methods):
		ls_i = lss[ni]
		cl_i = clrs[ni]

		wftype, diff_method = m

		print(f"----------------------------------")
		print(f">> Approx: {approx}")
		print(f">> diff_method: {diff_method}")
		print(f">> diff_order: {diff_order}")
		print(f"----------------------------------")

		Diff_values, time_diff = gw.get_derivatives(FreeParams,approx,GwPrms,dets,freq, jitgrad=False,
													diff_order=diff_order,diff_method=diff_method,
													step_size=1.e-5, enable_jax_waveforms=wftype=="jax",
													EarthRotation=True)
		time_vec, dt_total= time_diff

		print("Diff_Values:\n",np.shape(Diff_values))
		print("\nRuntimes:")

		'''
		for n, t in enumerate(time_vec):
			print(f"\t {FreeParams[n]}:\t dt = {t:.2e} seconds")
		print(f"\ntotal_time: {dt_total:.2e} seconds")
		'''

		print("shape(Diff_values):",np.shape(Diff_values))


		try:
			Diff_values = Diff_values.reshape(-1, 1000)
		except:
			Diff_values = np.array(Diff_values)

		ndim = int( np.ceil( np.sqrt(len(Diff_values)) ) )

		print('\nplotting_diffs...')
		N = len(Diff_values)

		Abs = lambda x: np.real( np.sqrt(x*np.conj(x)) )

		if diff_order == "first": 		labels = FreeParams
		elif diff_order == "second": 	labels = [ f"{x1}-{x2}" for x1 in FreeParams for x2 in FreeParams ]
		elif diff_order == "third": 	labels = [ f"{x1}-{x2}-{x3}" for x1 in FreeParams for x2 in FreeParams for x3 in FreeParams ]

		for i in trange(len(Diff_values)):
			plt.subplot(3,4,i+1)
			#plt.subplot(2,1,i+1)
			#plt.subplot(6,5,i+1)
			y = Abs(Diff_values[i])
			if diff_method != "autodiff":
				plt.loglog(freq,y,ls=ls_i,color=cl_i,lw=2,zorder=1)
			else:
				plt.loglog(freq,y,ls=ls_i,color=cl_i,lw=3,zorder=2)

			plt.ylim(.01*max(y),2*max(y))
			plt.title(labels[i])
			plt.xticks([])
			plt.yticks([])
			#plt.legend(fontsize=8)
			#plt.grid(which='both',alpha=0.3)

plt.plot([],'k-',label="jax-autodiff")

plt.suptitle(approx,weight="bold")
plt.tight_layout()
fig.savefig(f'outputs/fig_get_derivatives_{FreeParams[0]}-{FreeParams[1]}.png')
print("\nOK!")

plt.show()

print("\n",sys.argv[0],"OK!")