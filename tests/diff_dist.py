import sys, os
import GWDALI_v1 as gw
import numpy as np
import pandas as pd
import jax.numpy as jnp
import matplotlib.pyplot as plt

from time import time as now
from tqdm import trange

np.random.seed(0)

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

deltaM = np.sqrt(1-4*eta)
ln_Mc = np.log(Mc)
ln_eta = np.log(eta)
inv_eta = 1./eta

mass_keys = "Mc,eta,m1,m2,M,q,ln_Mc,ln_eta,inv_eta,deltaM".split(',')
values = Mc,eta,m1,m2,M,q,ln_Mc,ln_eta,inv_eta,deltaM

dL_prms = {}
dL_prms["dL"] = dL
dL_prms["inv_dL"] = 1./dL
dL_prms["inv_dL2"] = 1./dL**2
dL_prms["inv_sqrtdL"] = 1./np.sqrt(dL)
dL_prms["inv_lnDL"] = 1./np.log(dL)
dL_prms["lnDL"] = np.log(dL)

GwPrms = {}
GwPrms["ln_Mc"] = ln_Mc
GwPrms["deltaM"] = deltaM
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

approx = "TaylorF2"
diff_method = "numdiff"
diff_order  = "first"

print(f">> Running {approx}...\n")

fig = plt.figure(figsize=(8,6))
plt.rcParams.update({"font.family":"serif","font.size":20,"font.weight":"bold"})


for i, key in enumerate( dL_prms.keys() ):
	plt.subplot(2,3,i+1)
	FreeParams = [key]

	GwPrmsX = GwPrms.copy()
	GwPrmsX[key] = dL_prms[key]

	#print(GwPrmsX.keys()) ; quit()

	dets = [detector, detector]

	wftype, diff_method = "lal", "numdiff"
	#wftype, diff_method = "jax", "autodiff"

	#aux = "_EarthRotation"
	aux = ""

	Diff_values, time_diff = gw.get_derivatives(FreeParams,approx+aux,GwPrmsX,dets,freq, jitgrad=False,
												diff_order=diff_order,diff_method=diff_method,
												step_size=1.e-5, enable_jax_waveforms=wftype=="jax")
	time_vec, dt_total= time_diff

	#print("Diff_Values:\n",np.shape(Diff_values))
	#print("\nRuntimes:")
	#print("shape(Diff_values):",np.shape(Diff_values))

	try:	Diff_values = Diff_values.reshape(-1, 1000)
	except:	Diff_values = np.array(Diff_values)

	ndim = int( np.ceil( np.sqrt(len(Diff_values)) ) )

	print('\nplotting_diffs...')
	N = len(Diff_values)

	Abs = lambda x: np.real( np.sqrt(x*np.conj(x)) )

	if diff_order == "first":
		labels = FreeParams

	ls = {"jax":'-',"lal":"--"}
	cl = {"jax":'k',"lal":"r"}

	y = Abs(Diff_values[0])
	if diff_method != "autodiff":
		plt.loglog(freq,y,ls=ls[wftype],color=cl[wftype],lw=2)
	else:
		plt.loglog(freq,y,ls='-',color='g',lw=1)
	#plt.ylim(min(y[y!=0])*100,max(y)*10)
	plt.title(labels[0])
	plt.xticks([])
	plt.yticks([])

plt.suptitle(approx,weight="bold")
plt.plot([],[],'k-',label='jax [numdiff]')
plt.plot([],[],'r-',label='lal [numdiff]')
plt.plot([],[],'g-',label='jax [autodiff]')
fig.legend(loc='lower right')
plt.tight_layout()
fig.savefig(f'outputs/diffs/fig_get_derivatives_{FreeParams[0]}.png')
print("\nOK!")

plt.show()

print("\n",sys.argv[0],"OK!")