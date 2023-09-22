=================================
Usage
=================================

Example:

.. code:: python

    import numpy as np
    #-------------------
    import GWDALI as gw
    #-------------------
    from tqdm import trange
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(70,0.3)

    rad = np.pi/180 ; deg = 1./rad
    #--------------------------------------------
    # Detector, position and orientation
    #--------------------------------------------
    FreeParams = ['DL','iota','psi','phi_coal']

    # Cosmic Explorer:
    det0 = {"name":"CE","lon":-119,"lat":46,"rot":45,"shape":90}
    # Einstein Telescope:
    det1 = {"name":"ET","lon":10,"lat":43,"rot":0,"shape":60}
    det2 = {"name":"ET","lon":10,"lat":43,"rot":120,"shape":60}
    det3 = {"name":"ET","lon":10,"lat":43,"rot":-120,"shape":60}

    #------------------------------------------------------
    # Setting Injections (Single detection)
    #------------------------------------------------------
    z = 0.1 # Redshift

    params = {}
    params['m1']  = 1.3*(1+z) # mass of the first object [solar mass]
    params['m2']  = 1.5*(1+z) # mass of the second object [solar mass]
    params['z']   = z
    params['RA']       = np.random.uniform(-180,180)
    params['Dec']      = (np.pi/2-np.arccos(np.random.uniform(-1,1)))*deg
    params['DL']       = cosmo.luminosity_distance(z).value/1.e3 # Gpc
    params['iota']     = np.random.uniform(0,np.pi)      # Inclination angle (rad)
    params['psi']      = np.random.uniform(0,np.pi) # Polarization angle (rad)
    params['t_coal']   = 0  # Coalescence time
    params['phi_coal'] = 0  # Coalescence phase
    # Spins:
    params['sx1'] = 0 
    params['sy1'] = 0
    params['sz1'] = 0
    params['sx2'] = 0
    params['sy2'] = 0
    params['sz2'] = 0

    #----------------------------------------------------------------------
    # "approximant" options: 
    #       [Leading_Order, TaylorF2_py, ...] or any lal approximant
    #----------------------------------------------------------------------
    # "dali_method" options:
    #       [Fisher, Fisher_Sampling, Doublet, Triplet, Standard]
    #----------------------------------------------------------------------
    res = gw.GWDALI( Detection_Dict = params, 
             FreeParams     = FreeParams, 
             detectors      = [det0,det1,det2,det3], # Einstein Telescope + Cosmic Explorer
             approximant    = 'TaylorF2_py',
             dali_method    = 'Doublet',
             sampler_method = 'nestle', # Same as Bilby sampling method
             save_fisher    = False,
             save_cov       = False,
             plot_corner    = False,
             save_samples   = False,
             hide_info      = True,
             index          = 1,
             rcond          = 1.e-4,
             npoints=300) # points for "nested sampling" or steps/walkers for "MCMC"

    Samples = res['Samples']
    Fisher  = res['Fisher']
    CovFish = res['CovFisher']
    Cov     = res['Covariance']
    Rec = res['Recovery']
    Err     = res['Error']
    SNR     = res['SNR']