import numpy as np
import matplotlib.pyplot as plt
from time import time as now
from datetime import datetime
from getdist import MCSamples, plots

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
wf_type = "jax"

ndim = len(FreeParams)
sampler = "nestle"
npoints = 300
nsamples = 6000
ntemps = 6
nwalkers = 12
#============#============#============#============

tns_data = np.load("tns_outputs/tensors_lal_TaylorF2_numdiff.npz")
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
#============#============#============#============

os.system('rm -R outputs_bilby/')

samples_vec = []
for method in ["Fisher","Singlet","FisherG"]:
    if method=="FisherG":
        method = "Fisher"
        diagG = [ 1./10**2, 1./np.pi**2 ]
        FisherPrior = np.eye(2)*diagG
        dali_tensors = [Fisher+FisherPrior,Doublet,Triplet]

    t1 = now() ; print("Running GWDALI")
    res = gw.GWDALI(	GwPrms=GwPrms,
                        detectors=detectors,
                        FreeParams=FreeParams,
                        approx = approx,
                        method=method,
                        sampler=sampler,
                        new_priors = None,
                        #diff_method = "numdiff", # 'numdiff' or 'autodiff'
                        dali_tensors = dali_tensors,
                        step_size = [1.e-4,1.e-3,1.e-2],
                        plot_signal=False,
                        npoints=npoints,
                        nwalkers=nwalkers,
                        pos0 = None,
                        npool = 1,
                        ntemps = ntemps,
                        remove_out=False,
                        verbose = True,
                        hide_info = False,
                        output_name=f'outputs/gwdali_output/',
                        nsamples=nsamples,
                        enable_jax_waveforms=wf_type=="jax",
                        #thin_by_nact=1,
                        #burn_in_nact=0,
                        #burn_in_fixed_discard=0,
                        )
    print("\n\n\t\tMCMC Concluded!!!\n\n")
    
    if method == "Singlet": samples, priors, likelihood_times, Evidence = res[0]
    elif method == "Fisher": samples, Cov, eps = res[0]
    
    samples_vec.append(samples)

dt = int(now()-t1)

samplesF, samplesS, samplesG = samples_vec

samples1 = MCSamples(samples=samplesF, label="Fisher", names=FreeParams)
samples2 = MCSamples(samples=samplesG, label="Fisher + G.Prior", names=FreeParams)
samples3 = MCSamples(samples=samplesS, label="Singlet", names=FreeParams)

g = plots.get_subplot_plotter()
g.triangle_plot([samples1, samples2,samples3], filled=True)

g.export(f"outputs/fig_Fisher-Singlet.jpg")
plt.show()

print("Fisher/Singlet OK!")