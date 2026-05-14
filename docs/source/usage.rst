=================================
Usage
=================================

Example:

.. code:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from time import time as now
    from datetime import datetime
    from corner import corner
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
    GwPrms['dL']    = 10.96
    GwPrms['Mc']    = 9.57
    GwPrms['eta']   = 0.25
    GwPrms['iota']  = 0.52
    GwPrms['psi']   = 2.32
    GwPrms['t_coal']    = 0.00
    GwPrms['phi_coal']  = -1.37
    GwPrms['RA']    = -135.75
    GwPrms['Dec']   = 72.92
    GwPrms['sz1']   = 0.31
    GwPrms['sz2']   = 0.16

    GwPrms['sx1']   = 0.00
    GwPrms['sy1']   = 0.00

    GwPrms['sx2']   = 0.00
    GwPrms['sy2']   = 0.00

    #============#============#============#============
    FreeParams = f"dL,iota".split(',')
    ndim = len(FreeParams)

    approx = "TaylorF2"
    method = "Doublet"
    wf_type = "jax"

    ndim = len(FreeParams)
    sampler = "nestle"
    npoints = 300
    nsamples = 6000
    ntemps = 6
    nwalkers = 12

    os.system('rm -R outputs_bilby/')
    t1 = now() ; print("Running GWDALI")
    res = gw.GWDALI(    GwPrms=GwPrms,
                        detectors=detectors,
                        FreeParams=FreeParams,
                        approx = approx,
                        method=method,
                        sampler=sampler,
                        new_priors = None,
                        diff_method = "numdiff", # 'numdiff' or 'autodiff'
                        # dali_tensors = dali_tensors, 
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
                        output_name=f'gwdali_output/',
                        nsamples=nsamples,
                        enable_jax_waveforms=wf_type=="jax",
                        #thin_by_nact=1,
                        #burn_in_nact=0,
                        #burn_in_fixed_discard=0,
                        )
    print("\n\n\t\tMCMC Concluded!!!\n\n")
    Results, Truths, Tensors, Fisher_Matrix, Time = res
    time_dali, time_mcmc = Time

    samples, priors, likelihood_times, Evidence = Results

    logZ, logZ_err = Evidence

    dt = int(now()-t1)

    file_name = f"samples_{wf_type}_{method}.txt"
    np.savetxt(file_name,samples,fmt='%e',delimiter='\t',header=f'Concluded in {dt} seconds\ndL, iota')

    fig = corner(samples,bins=40,color='k',smooth=1,smooth1d=1,fill_contours=True,labels=FreeParams)
    plt.plot([],[],'k-',label=f"{wf_type}_{method}")
    fig.legend(loc='upper right')
    plt.suptitle(approx,weight="bold",ha="right")
    fig.savefig(f"fig_{wf_type}_{method}.jpg")
    plt.show()