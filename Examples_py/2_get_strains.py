import GWDALI as gw
import numpy as np
import matplotlib.pyplot as plt
from time import time as now

Approxs = ["TaylorF2"] + [f"IMRPhenom{x}" for x in "A,B,C,D,HM".split(',')]

rad = np.pi/180
deg = 1./rad

detector = {'name':'ET','lon':0, 'lat':90, 'rot':0, 'shape':90}
detectors = [detector]

GwPrms = {}
GwPrms['dL']       = np.random.uniform(1,10) # Gpc
GwPrms['Mc']       = np.random.uniform(1,100)
GwPrms['eta']      = np.random.uniform(0.08,0.25)
GwPrms['iota']     = np.random.uniform(0,np.pi)
GwPrms['psi']      = np.random.uniform(0,np.pi)
GwPrms['t_coal']   = np.random.uniform(-1,1)
GwPrms['phi_coal'] = np.random.uniform(0.,np.pi)
GwPrms['RA']       = np.random.uniform(0,360)
GwPrms['Dec']      = np.random.uniform(-90,90)

GwPrms["sx1"]  = 0.
GwPrms["sy1"]  = 0.
GwPrms["sz1"]  = np.random.uniform(-0.98,0.98)

GwPrms["sx2"]  = 0.
GwPrms["sy2"]  = 0.
GwPrms["sz2"]  = np.random.uniform(-0.98,0.98)

sz1 = GwPrms['sz1']
sz2 = GwPrms['sz2']
eta = GwPrms['eta']

dM = np.sqrt(1.-4*eta)
chi_s = .5*(sz1+sz2)
chi_a = .5*(sz1-sz2)
chi_eff = chi_s + dM*chi_a

fig = plt.figure(figsize=(12,10))
plt.suptitle("GWDALI 1.0: Strain ($h=F_+h_++F_{\\times}h_{\\times}$) [JAX vs LAL]",weight="bold")

Amplitude = lambda x: np.real( np.sqrt(x*np.conj(x)) )

Diff = lambda x,y: np.abs((x-y)/y)

for n, approx in enumerate(Approxs):
	print(f">> Running {approx}...\n")
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

	
	H_vec_jax = gw.get_strain(detectors,GwPrms,freq,approx,enable_jax_waveforms=True,disable_jit=True,EarthRotation=False)
	H_vec_lal = gw.get_strain(detectors,GwPrms,freq,approx,enable_jax_waveforms=False,dF=dF,EarthRotation=False)

	hj = H_vec_jax[0] # first detector 
	hl = H_vec_lal[0]

	ax = plt.subplot(6,1,n+1)
	plt.xticks([]) ; plt.yticks([])
	for spine in ax.spines.values():
		spine.set_visible(False)
		
	plt.subplot(6,2,2*n+1)

	AmpJ = Amplitude(hj)
	AmpL = Amplitude(hl)

	PhaseJ = np.abs( np.unwrap(np.angle(hj/AmpJ)) ) / np.pi
	PhaseL = np.abs( np.unwrap(np.angle(hl/AmpL)) ) / np.pi

	plt.plot([],'w.',label=approx)
	plt.loglog(freq,AmpJ,color='k',lw=1,ls='-')
	plt.loglog(freq,AmpL,color='r',lw=3,ls=':')
	plt.ylim(1.e-26,1.e-21)
	plt.grid(ls='--',which='both',alpha=0.3)
	plt.ylabel("Amplitude")
	plt.legend(loc='upper right')
	if (n!= len(Approxs)-1): plt.xticks([])
	else: plt.xlabel('frequency [Hz]')

	plt.subplot(6,2,2*n+2)
	plt.ylabel(f"| Phase [$\\pi$] |")
	plt.plot(freq,PhaseJ,color='k',lw=1,ls='-')
	plt.plot(freq,PhaseL,color='r',lw=3,ls=':')
	plt.xscale('log') ; plt.yscale('log')
	plt.grid(ls='--',which='both',alpha=0.3)
	if (n!= len(Approxs)-1): plt.xticks([])
	else: plt.xlabel('frequency [Hz]')

	print(approx,"OK!")

plt.tight_layout()
plt.subplots_adjust(right=.85)
plt.plot([],'k-',label='jax')
plt.plot([],'r:',label='lal')
for key in GwPrms.keys():
	plt.plot([],[],'k.',label=f"{key} = {GwPrms[key]:.2f}")
fig.legend(loc='upper right')
fig.savefig('outputs/fig_get_strains.png')

plt.show()