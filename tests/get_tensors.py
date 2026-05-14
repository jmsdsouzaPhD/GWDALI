import sys, os
import GWDALI_v1 as gw
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

dL = 450.

det_ET1 = {'name':'ET','lon':6, 'lat':50, 'rot':0, 'shape':60}
det_ET2 = {'name':'ET','lon':6, 'lat':50, 'rot':120, 'shape':60}
det_ET3 = {'name':'ET','lon':6, 'lat':50, 'rot':-120, 'shape':60}
detectors = [det_ET1, det_ET2, det_ET3]

freq = 10**np.linspace(np.log10(5),3,1000)

DataFrameSrc = pd.read_csv(f"300_srcs.csv")

#=========================================================================================
Idx = 0
dL = DataFrameSrc["dL"][Idx]

GwPrms = {}
GwPrms['dL'] = dL
GwPrms['Mc']   = DataFrameSrc["Mc"][Idx]
GwPrms['eta']  = DataFrameSrc["eta"][Idx]
GwPrms['iota'] = DataFrameSrc["iota"][Idx]
GwPrms['psi']  = DataFrameSrc["psi"][Idx]
GwPrms['t0']   = DataFrameSrc["t0"][Idx]
GwPrms['phi0'] = DataFrameSrc["phi0"][Idx]
GwPrms['RA']   = DataFrameSrc["RA"][Idx]
GwPrms['Dec']  = DataFrameSrc["Dec"][Idx]
GwPrms["sz1"]  = DataFrameSrc["sz1"][Idx]
GwPrms["sz2"]  = DataFrameSrc["sz2"][Idx]

GwPrms["sx1"] = 0.
GwPrms["sy1"] = 0.
GwPrms["sx2"] = 0.
GwPrms["sy2"] = 0.

fig = plt.figure(figsize=(12,12))

approx = "TaylorF2"#_Spins_ISCO"

print(f">> Running {approx}...\n")

FreeParams = f"dL,iota".split(',')

res = gw.get_dali_tensors(	GwPrms,
							detectors,
							FreeParams=["dL","iota"],
							method="Doublet",
							approx="TaylorF2",
							enable_jax_waveforms=False,
							save_tensors=["tns_outputs","src1"],
						    diff_method="numdiff",
						    step_size=[1.e-6],
						    hide_info=False)

Tensors, time_dets_tns = res

Fisher = Tensors["Fisher"]
Db12, Db22 = Tensors["Doublet"]
Tp13, Tp23, Tp33 = Tensors["Triplet"]

#np.savez_compressed(f"outputs/tensors_{wf_type}_{approx}_{diff_method}.npz",Fisher=Fisher,Db12=Db12,Db22=Db22,Tp13=Tp13,Tp23=Tp23,Tp33=Tp33)
print("Tensors Saved!!!")
