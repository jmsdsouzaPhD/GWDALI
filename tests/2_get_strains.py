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

freq = 10**np.linspace(np.log10(5),3,10000)

Mc = 44.#36.
eta = 0.24
dL = 460.
iota = np.pi/3
m1 = 0.5*(Mc/eta**(3./5)) * (1 + np.sqrt(1-4*eta))
m2 = 0.5*(Mc/eta**(3./5)) * (1 - np.sqrt(1-4*eta))
M = m1+m2
sz1 = 0.5
sz2 = -0.2

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

fig = plt.figure(figsize=(12,10))
plt.suptitle("GWDALI 1.0: Strain ($h=F_+h_++F_{\\times}h_{\\times}$) [JAX vs LAL]",weight="bold")

Amplitude = lambda x: np.real( np.sqrt(x*np.conj(x)) )
Phase = lambda x: np.unwrap(np.angle(x))
Diff = lambda x,y: np.abs((x-y)/y)

for n, approx in enumerate(Approxs):
	print(f">> Running {approx}...\n")
	if approx == "IMRPhenomA": GwPrms["sz1"] = GwPrms["sz2"]=sz1=sz2=0

	chi1, chi2 = GwPrms["sz1"], GwPrms["sz2"]
	chi_eff = (chi1*m1 + chi2*m2)/M
	iota = GwPrms["iota"]
	fRef = 1.0
	dF = 0.01

	for ni, enable_jax_waveforms in enumerate([True,False]):
		t1 = now()
		if approx == "IMRPhenomC": GwPrms["sz1"] = GwPrms["sz2"] = chi1 = chi2 = chi_eff
		if enable_jax_waveforms:
			H_vec = gw.get_strain(detectors,GwPrms,freq,approx,enable_jax_waveforms,disable_jit=True,EarthRotation=True)
		else:
			H_vec = gw.get_strain(detectors,GwPrms,freq,approx,enable_jax_waveforms,dF=dF,EarthRotation=True)

		h, hr, hi = [x[0] for x in H_vec]

		dt = now()-t1

		color = ['k','r'][ni]
		ls = ['-',':'][ni]
		lw = [1,2][ni]
		label = ['jax','lal'][ni]

		ax = plt.subplot(6,1,n+1)
		plt.title(approx,weight="bold")
		plt.xticks([]) ; plt.yticks([])
		for spine in ax.spines.values():
			spine.set_visible(False)
			
		plt.subplot(6,2,2*n+1)

		plt.loglog(freq,Amplitude(h),color=color,lw=lw,ls=ls,label=label)
		plt.ylim(1.e-26,1.e-21)
		plt.grid(ls='--',which='both',alpha=0.3)
		#plt.ylim(h_amp[0]/1000,h_amp[0]*2)
		plt.legend(loc='lower left')
		plt.ylabel("Amplitude")
		if (n!= len(Approxs)-1): plt.xticks([])
		else: plt.xlabel('frequency [Hz]')

		plt.subplot(6,2,2*n+2)
		plt.ylabel(f"Phase")
		plt.plot(freq,Phase(h),color=color,lw=lw,ls=ls)
		plt.xscale('log')
		plt.grid(ls='--',which='both',alpha=0.3)
		#plt.ylim(h_amp[0]/1000,h_amp[0]*2)
		#plt.legend(loc='upper right')
		if (n!= len(Approxs)-1): plt.xticks([])
		else: plt.xlabel('frequency [Hz]')

		print(approx,f"[enable_jax_waveforms: {enable_jax_waveforms}] OK!")

plt.tight_layout()
fig.savefig('outputs/fig_get_strains.png')

plt.show()

print("\n",sys.argv[0],"OK!")