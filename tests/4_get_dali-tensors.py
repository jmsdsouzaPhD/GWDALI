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

dL = 450.

det_ET1 = {'name':'ET','lon':6, 'lat':50, 'rot':0, 'shape':60}
det_ET2 = {'name':'ET','lon':6, 'lat':50, 'rot':120, 'shape':60}
det_ET3 = {'name':'ET','lon':6, 'lat':50, 'rot':-120, 'shape':60}
detectors = [det_ET1, det_ET2, det_ET3]

freq = np.logspace(1,3,10000)

GwPrms = {}                                                                                                
GwPrms['dL'] 	= 10.96
GwPrms['Mc'] 	= 9.57
GwPrms['eta'] 	= 0.25
GwPrms['iota'] 	= 0.52
GwPrms['psi'] 	= 2.32
GwPrms['t_coal'] 	= 0.00
GwPrms['phi_coal'] 	= -1.37
GwPrms['RA'] 	= -135.75
GwPrms['Dec'] 	= 72.92
GwPrms['sz1'] 	= 0.31
GwPrms['sz2'] 	= 0.16

GwPrms['sx1'] 	= 0.00
GwPrms['sy1'] 	= 0.00

GwPrms['sx2'] 	= 0.00
GwPrms['sy2'] 	= 0.00

approx = "TaylorF2"
if approx == "IMRPhenomA": GwPrms["sz1"] = GwPrms["sz2"] = 0. 

print(f">> Running {approx}...\n")

FreeParams = f"dL,iota".split(',')

for wf_type in ["lal","jax"]:
	if wf_type == "lal":
		DiffMethods = ["numdiff"]
	else:
		DiffMethods = ["numdiff"]#,"autodiff"]

	for diff_method in DiffMethods:
		res = gw.get_dali_tensors(GwPrms,detectors,FreeParams,"Triplet",approx,
								  enable_jax_waveforms=wf_type=="jax",
								  #save_tensors=["tns_outputs","src1"],
							      diff_method=diff_method,
							      #step_size=[],
							      hide_info=False)

		Tensors, time_dets_tns = res

		Fisher = Tensors["Fisher"]
		Db12, Db22 = Tensors["Doublet"]
		Tp13, Tp23, Tp33 = Tensors["Triplet"]

		np.savez_compressed(f"tns_outputs/tensors_{wf_type}_{approx}_{diff_method}.npz",Fisher=Fisher,Db12=Db12,Db22=Db22,Tp13=Tp13,Tp23=Tp23,Tp33=Tp33)
		print("Tensors Saved!!!",wf_type,diff_method)
