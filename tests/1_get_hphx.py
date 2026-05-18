import sys, os
import GWDALI as gw
import numpy as np
import pandas as pd
import jax.numpy as jnp
import matplotlib.pyplot as plt

from astropy.cosmology import FlatLambdaCDM
from time import time as now
from tqdm import trange

Approxs = ["TaylorF2"] + [f"IMRPhenom{x}" for x in "A,B,C,D,HM".split(',')]

N = len(Approxs)

rad = np.pi/180
deg = 1./rad

GwPrms = {}
GwPrms['dL']       = np.random.uniform(1,10) # Gpc
GwPrms['Mc']       = np.random.uniform(1,100)
GwPrms['eta']      = np.random.uniform(0.08,0.25)
GwPrms['iota']     = np.random.uniform(0,np.pi)
GwPrms['t_coal']   = np.random.uniform(-1,1)
GwPrms['phi_coal'] = np.random.uniform(0.,np.pi)

GwPrms["sx1"]  = 0.
GwPrms["sy1"]  = 0.
GwPrms["sz1"]  = np.random.uniform(-0.98,0.98)

GwPrms["sx2"]  = 0.
GwPrms["sy2"]  = 0.
GwPrms["sz2"]  = np.random.uniform(-0.98,0.98)

Mc  = GwPrms["Mc"]
eta = GwPrms["eta"]
sz1 = GwPrms["sz1"]
sz2 = GwPrms["sz2"]

dM = np.sqrt(1-4*eta)
M = Mc/eta**(3./5)
m1 = .5*M*(1+dM)
m2 = .5*M*(1-dM)
chi_eff = (sz1*m1 + sz2*m2)/M

fig = plt.figure(figsize=(16,N*2))
plt.suptitle("GWDALI 1.0: GW Polarization [$h \\equiv h_+-ih_{\\times}$] (Jax vs LAL)",weight="bold")

Amplitude = lambda x: np.real( np.sqrt(x*np.conj(x)) )
def Phase(y):
	A = y/Amplitude(y)
	Cos, Sin = A.real, A.imag
	return np.arctan2(Sin,Cos) / np.pi

Diff = lambda x,y: np.abs((x-y)/y)

x = np.logspace(1,3,1000)
y = x**(-7./6)*np.exp(1.j*( 2*np.pi*x + 1.e2*x**(-5./3)) )

def NotNan(x):
	y = []
	for xi in x:
		if (not np.isnan(xi)) and (not np.isinf(xi)):
			y.append(xi)
	return np.array(y)

for n, approx in enumerate(Approxs):
	print(f">> Running {approx}...\n")
	if approx == "IMRPhenomA": 
		GwPrms["sz1"], GwPrms["sz2"] = 0. , 0.
	elif approx == "IMRPhenomC":
		GwPrms["sz1"], GwPrms["sz2"] = chi_eff, chi_eff
	else:
		GwPrms["sz1"], GwPrms["sz2"] = sz1, sz2
	
	fRef = 1.0
	dF = 0.01

	freq = np.arange(1.,1.e3,dF)
	disable_jit = True

	hphx_vec_lal = gw.get_hphx(GwPrms,freq,approx,enable_jax_waveforms=False,dF=dF)
	hphx_vec_jax = gw.get_hphx(GwPrms,freq,approx,enable_jax_waveforms=True,disable_jit=disable_jit,EarthRotation=False)

	hp_lal, hx_lal = hphx_vec_lal
	hp_jax, hx_jax = hphx_vec_jax

	h_jax = hp_jax - 1.j*hx_jax
	h_lal = hp_lal - 1.j*hx_lal

	amp_jax = Amplitude(h_jax)
	phi_jax = Phase(h_jax)

	amp_lal = Amplitude(h_lal)
	phi_lal = Phase(h_lal)

	diffA = Diff(amp_jax,amp_lal)
	diffP = np.abs(phi_jax-phi_lal) # units of pi

	ax = plt.subplot(N,1,n+1)
	plt.title(approx,weight="bold")
	plt.xticks([]) ; plt.yticks([])
	for spine in ax.spines.values():
		spine.set_visible(False)

	plt.subplot(N,4,4*n+1)
	plt.loglog(freq,amp_jax,'k-')
	plt.loglog(freq,amp_lal,'r:')
	plt.ylim(amp_jax[-1],amp_jax[0])
	plt.xscale('log')
	plt.grid(which='both',alpha=.3)
	plt.xlim(freq[1],freq[-1])
	plt.ylabel("Amplitude")
	if (n!= len(Approxs)-1): plt.xticks([])
	else: plt.xlabel('frequency [Hz]')

	plt.subplot(N,4,4*n+3)
	plt.loglog(freq[1:],diffA[1:],'b.',ms=1)
	plt.xscale('log')
	plt.ylabel("$\\Delta\\ Amp.$")
	plt.grid(which='both',alpha=.3)
	plt.xlim(freq[1],freq[-1])
	if (n!= len(Approxs)-1): plt.xticks([])
	else: plt.xlabel('frequency [Hz]')

	plt.subplot(N,4,4*n+2)
	plt.plot(freq,phi_jax,'k-')
	plt.plot(freq,phi_lal,'r:')
	plt.xlim(freq[1],freq[-1])
	plt.ylabel("Phase [$\\pi$]")

	plt.xscale('log')
	plt.grid(which='both',alpha=.3)
	if (n!= len(Approxs)-1): plt.xticks([])
	else: plt.xlabel('frequency [Hz]')

	plt.subplot(N,4,4*n+4)
	plt.loglog(freq,diffP,'b.',ms=1)
	plt.xscale('log')
	plt.xlim(freq[1],freq[-1])
	plt.grid(which='both',alpha=.3)
	plt.ylabel("$\\Delta\\ Phase$ [$\\pi$]")
	if (n!= len(Approxs)-1): plt.xticks([])
	else: plt.xlabel('frequency [Hz]')

	print(approx,"OK!")

plt.tight_layout()

plt.tight_layout()
plt.subplots_adjust(right=.85)
plt.plot([],'k-',label='jax')
plt.plot([],'r:',label='lal')
for key in GwPrms.keys():
	plt.plot([],[],'k.',label=f"{key} = {GwPrms[key]:.2f}")
fig.legend(loc='upper right')

fig.savefig('outputs/fig_get_hphx.pdf')

print("\n",sys.argv[0],"OK!")
plt.show()