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
FreeParams = ['DL','iota','psi']


rot_LHO = 180+36
rot_LLO = 270+18
rot_Virgo = 180-19
rot_KAGRA = 90+28.31

det1 = {"name":"aLIGO","lon":-119.4073,"lat":46.4554,"rot":rot_LHO,"shape":90} # Hanford
det2 = {"name":"aLIGO","lon":-90.7740,"lat":30.5629,"rot":rot_LLO,"shape":90} # Livingston
det3 = {"name":"aVirgo","lon":10.5024,"lat":43.6279,"rot":rot_Virgo,"shape":90}
det4 = {"name":"KAGRA","lon":137.1761,"lat":36.2539,"rot":rot_KAGRA,"shape":90}

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
params['iota']     = np.random.uniform(0,np.pi)          # Inclination angle (rad)
params['psi']      = np.random.uniform(-np.pi,np.pi) # Polarization angle (rad)
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
#               [Leading_Order, TaylorF2_py] or any lal approximant
#----------------------------------------------------------------------
# "dali_method" options:
#               [Fisher, Fisher_Sampling, Doublet, Triplet, Standard]
#----------------------------------------------------------------------
res = gw.GWDALI( Detection_Dict = params,
                 FreeParams     = FreeParams,
                 detectors      = [det1,det2,det3,det4], # LIGO-H, LIGO-L, Virgo and KAGRA
                 approximant    = 'TaylorF2_py',
                 dali_method    = 'Fisher_Sampling',
                 sampler_method = 'nestle', # Same as Bilby sampling method
                 save_fisher    = True,
                 save_cov       = True,
                 plot_corner    = True,
                 save_samples   = True,
                 hide_info      = True,
                 index          = 1,
                 rcond          = 1.e-4,
                 npoints=50) # points for "nested sampling" or steps/walkers for "MCMC"

Samples = res['Samples']
Fisher  = res['Fisher']
CovFish = res['CovFisher']
Cov     = res['Covariance']
Rec     = res['Recovery']
Err     = res['Error']
SNR     = res['SNR']
