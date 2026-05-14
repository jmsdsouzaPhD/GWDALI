import sys, os
import GWDALI as gw
import numpy as np
import pandas as pd
import jax.numpy as jnp
import matplotlib.pyplot as plt

from astropy.cosmology import FlatLambdaCDM
from time import time as now
from tqdm import trange

np.random.seed(0)

Approxs = ["TaylorF2"] + [f"IMRPhenom{x}" for x in "A,B,C,D,HM".split(',')]
N = len(Approxs)

rad = np.pi/180
deg = 1./rad

det_CE = {'name':'CE','lon':0, 'lat':90, 'rot':0, 'shape':90}
det_ET1 = {'name':'ET','lon':0, 'lat':90, 'rot':0, 'shape':60}
det_ET2 = {'name':'ET','lon':0, 'lat':90, 'rot':120, 'shape':60}
det_ET3 = {'name':'ET','lon':0, 'lat':90, 'rot':-120, 'shape':60}
detectors = [det_CE]#det_ET1, det_ET2, det_ET3]

#freq = np.logspace(1,4,10000)
freq = 10**np.linspace(np.log10(5),3,1000)

Mc = 50.#44.#36.
eta = 0.2

dL = 460.
iota = np.pi/3
m1 = 0.5*(Mc/eta**(3./5)) * (1 + np.sqrt(1-4*eta))
m2 = 0.5*(Mc/eta**(3./5)) * (1 - np.sqrt(1-4*eta))
M = m1+m2
sz1 = 0.5
sz2 = -0.2

#print(m1,m2) ; quit()

GwPrms = {}
GwPrms['dL']       = dL/1.e3 # Gpc
GwPrms['Mc']       = np.random.uniform(1,100)
GwPrms['eta']      = np.random.uniform(0.08,0.25)
GwPrms['iota']     = np.random.uniform(0,np.pi)
GwPrms['psi']      = np.random.uniform(0,np.pi)
GwPrms['t_coal']   = np.random.uniform(-1,1)
GwPrms['phi_coal'] = np.random.uniform(0.,np.pi)
GwPrms['RA']       = np.random.uniform(-180,180)
GwPrms['Dec']      = np.random.uniform(-90,90)

GwPrms["sx1"]  = 0.
GwPrms["sy1"]  = 0.
GwPrms["sz1"]  = 0*np.random.uniform(-0.98,0.98)

GwPrms["sx2"]  = 0.
GwPrms["sy2"]  = 0.
GwPrms["sz2"]  = 0*np.random.uniform(-0.98,0.98)

fig = plt.figure(figsize=(16,N*2))
plt.suptitle("GWDALI 1.0: GW Polarization [$h \\equiv h_+-ih_{\\times}$] (Jax vs LAL)",weight="bold")

Amplitude = lambda x: np.real( np.sqrt(x*np.conj(x)) )
Phase = lambda x: np.unwrap(np.angle(x))
Diff = lambda x,y: np.abs((x-y)/y)

for n, approx in enumerate(Approxs):
	print(f">> Running {approx}...\n")
	if approx == "IMRPhenomA": GwPrms["sz1"] = GwPrms["sz2"]=sz1=sz2=0

	chi1, chi2 = GwPrms["sz1"], GwPrms["sz2"]
	chi_eff = (chi1*m1 + chi2*m2)/M
	fRef = 1.0
	
	dF = 0.01
	disable_jit = True

	if approx == "IMRPhenomC": GwPrms["sz1"] = GwPrms["sz2"] = chi_eff

	hphx_vec_lal = gw.get_hphx(detectors,GwPrms,freq,approx,enable_jax_waveforms=False,dF=dF)
	hphx_vec_jax = gw.get_hphx(detectors,GwPrms,freq,approx,enable_jax_waveforms=True,disable_jit=disable_jit,EarthRotation=True)
	
	hp_lal, hx_lal = hphx_vec_lal[0]
	hp_jax, hx_jax = hphx_vec_jax[0]

	h_jax = hp_jax - 1.j*hx_jax
	h_lal = hp_lal - 1.j*hx_lal

	amp_jax = Amplitude(h_jax)
	phi_jax = Phase(h_jax)

	amp_lal = Amplitude(h_lal)
	phi_lal = Phase(h_lal)

	diffA = Diff(amp_jax,amp_lal)
	diffP = Diff(phi_jax,phi_lal)

	ax = plt.subplot(N,1,n+1)
	plt.title(approx,weight="bold")
	plt.xticks([]) ; plt.yticks([])
	for spine in ax.spines.values():
		spine.set_visible(False)

	plt.subplot(N,4,4*n+1)
	plt.loglog(freq,amp_jax,'k-')
	plt.loglog(freq,amp_lal,'r:',lw=3)
	plt.ylim(amp_jax[-1],amp_jax[0])
	#plt.plot(freq,h_lal,'k-')
	#plt.plot(freq,h_jax,'r:')
	plt.xscale('log')
	plt.grid(which='both',alpha=.3)
	plt.xlim(freq[1],freq[-1])
	plt.ylabel("Amplitude")
	if (n!= len(Approxs)-1): plt.xticks([])
	else: plt.xlabel('frequency [Hz]')

	plt.subplot(N,4,4*n+3)
	plt.loglog(freq[1:],diffA[1:],'b.',ms=1)
	#plt.plot(freq,h_lal,'k-')
	#plt.plot(freq,h_jax,'r:')
	plt.xscale('log')
	plt.ylabel("$\\Delta\\ Amp.$")
	plt.grid(which='both',alpha=.3)
	plt.xlim(freq[1],freq[-1])
	if (n!= len(Approxs)-1): plt.xticks([])
	else: plt.xlabel('frequency [Hz]')

	plt.subplot(N,4,4*n+2)
	plt.plot(freq,phi_jax,'k-')
	plt.plot(freq,phi_lal,'r:',lw=3)
	plt.xlim(freq[1],freq[-1])
	plt.ylabel("Phase")
	#plt.plot(freq,hx_lal,'k-')
	#plt.plot(freq,hx_jax,'r:')

	plt.xscale('log')
	plt.grid(which='both',alpha=.3)
	if (n!= len(Approxs)-1): plt.xticks([])
	else: plt.xlabel('frequency [Hz]')

	plt.subplot(N,4,4*n+4)
	plt.loglog(freq,diffP,'b.',ms=1)
	#plt.plot(freq,h_lal,'k-')
	#plt.plot(freq,h_jax,'r:')
	plt.xscale('log')
	plt.xlim(freq[1],freq[-1])
	plt.grid(which='both',alpha=.3)
	plt.ylabel("$\\Delta\\ Phase$")
	if (n!= len(Approxs)-1): plt.xticks([])
	else: plt.xlabel('frequency [Hz]')

	print(approx,"OK!")

plt.plot([],'k-',label='jax')
plt.plot([],'r:',label='lal')
fig.legend(loc='upper right')

plt.tight_layout()
fig.savefig('outputs/fig_get_hphx.png')
fig.savefig('outputs/fig_get_hphx.pdf')

print("\n","---"*10,'\n')
print("dF =", dF)
print("disable_jit =", disable_jit)
print("\n",sys.argv[0],"OK!")

plt.show()