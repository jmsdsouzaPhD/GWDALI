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
dL = 460. # Mpc
iota = np.pi/3
m1 = 0.5*(Mc/eta**(3./5)) * (1 + np.sqrt(1-4*eta))
m2 = 0.5*(Mc/eta**(3./5)) * (1 - np.sqrt(1-4*eta))
M = m1+m2
sz1 = 0.5
sz2 = -0.2

GwPrms = {}
GwPrms['dL']   = dL/1.e3 # Gpc
GwPrms['Mc']   = Mc
GwPrms['eta']  = eta
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


Abs = lambda x: np.real( np.sqrt(x*np.conj(x)) )

for n, approx in enumerate(Approxs):
	print(f"\n>> Running {approx}...")
	if approx == "IMRPhenomA": GwPrms["sz1"] = GwPrms["sz2"]=sz1=sz2=0

	chi1, chi2 = GwPrms["sz1"], GwPrms["sz2"]
	chi_eff = (chi1*m1 + chi2*m2)/M
	phi0, t0 = GwPrms["phi_coal"], GwPrms["t_coal"]
	iota = GwPrms["iota"]
	fRef = 1.0

	for ni, wf_type in enumerate(["jax","lal"]):
		enable_jax_waveforms = wf_type=="jax"
		t1 = now()
		if approx == "IMRPhenomC": GwPrms["sz1"] = GwPrms["sz2"] = chi1 = chi2 = chi_eff
		SNR_vec, SNR_tot = gw.get_SNR(detectors,GwPrms,approx,enable_jax_waveforms)
		
		print(f"\t[{wf_type}] SNR = {SNR_tot}",end='\t')
	print()