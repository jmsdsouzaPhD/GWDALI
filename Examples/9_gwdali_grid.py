import numpy as np
import matplotlib.pyplot as plt
from time import time as now
from datetime import datetime
from jcorner import jcorner
import traceback
import pandas as pd
import GWDALI as gw
import sys, os
import scipy

from tqdm import trange

det_ET1 = {'name':'ET','lon':6, 'lat':50, 'rot':0, 'shape':60}
det_ET2 = {'name':'ET','lon':6, 'lat':50, 'rot':120, 'shape':60}
det_ET3 = {'name':'ET','lon':6, 'lat':50, 'rot':-120, 'shape':60}
detectors = [det_ET1, det_ET2, det_ET3]

#=========================================================================================
rad = np.pi/180

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

#============#============#============#============
FreeParams = f"dL,iota".split(',')
ndim = len(FreeParams)

approx = "TaylorF2"
method = "Doublet"
#wf_type = "jax"
#wf_type = "lal"
wf_type = sys.argv[1] # jax or lal

ndim = len(FreeParams)
sampler = "grid"
npoints = 100

tns_data = np.load("../tns_outputs/tensors_lal_TaylorF2_numdiff.npz")
Fisher = tns_data["Fisher"]
Db12 = tns_data["Db12"]
Db22 = tns_data["Db22"]
Tp13 = tns_data["Tp13"]
Tp23 = tns_data["Tp23"]
Tp33 = tns_data["Tp33"]
Doublet = [Db12, Db22]
Triplet = [Tp13, Tp23, Tp33]
dali_tensors = [Fisher,Doublet,Triplet]
print("dali_tensors load successfully!")
#quit()
#============#============#============#============

limits = {}
limits['dL'] = np.linspace(1,3*GwPrms['dL'],npoints)
limits['iota'] = np.linspace(0*rad,100*rad,npoints)

t1 = now() ; print("Running GWJAX")
res = gw.GWDALI(	GwPrms=GwPrms,
					detectors=detectors,
					FreeParams=FreeParams,
					approx = approx,
					method=method,
					sampler=sampler,
					new_priors = None,
					diff_method = "numdiff", # 'numdiff' or 'autodiff'
					limits = limits,
					#dali_tensors = dali_tensors,
					step_size = [1.e-4,1.e-3,1.e-2],
					plot_signal=False,
					npoints=npoints,
					#nwalkers=nwalkers,
					pos0 = None,
					npool = 1,
					#ntemps = ntemps,
					remove_out=False,
					verbose = True,
					hide_info = False,
					output_name=f'gwdali_output/',
					#nsamples=nsamples,
					enable_jax_waveforms=wf_type=="jax",
					#thin_by_nact=1,
					#burn_in_nact=0,
					#burn_in_fixed_discard=0,
					)

Results, Truths, Tensors, Fisher_Matrix, Time = res
time_dali, time_mcmc = Time

os.system("clear")
Prms, Posterior, likelihood_times = Results
print(np.min(Posterior),np.max(Posterior))

xx, yy = np.meshgrid(Prms['dL'],Prms['iota'])

fig = plt.figure(figsize=(10,5))
plt.suptitle(f"{method} - {wf_type}",weight="bold")

plt.subplot(121)
c = plt.contourf(xx,yy,Posterior,cmap='magma',levels=30)
plt.plot(GwPrms['dL'],GwPrms['iota'],'wo',mec='k',label='Injection')
plt.xlabel('$d_L$ [Gpc]')
plt.ylabel('$\\iota$ [rad]')
plt.grid(alpha=.3)
plt.legend()
plt.colorbar(c,label='posterior',orientation='horizontal')

plt.subplot(122)
plt.hist(likelihood_times[1:]/1.e-3,bins=100,histtype='step',lw=2,color='k',density=True)
plt.xlabel('likelihood call time [ms]')
plt.grid()

plt.tight_layout()
#fig.savefig(f"outputs/fig_grid-samples_{wf_type}_{method}.png")
fig.savefig(f"outputs/fig_grid-samples_{wf_type}_{method}_no-tensors.png")

#fig.savefig(f"Grid_{wf_type}_{method}.png")

plt.show()