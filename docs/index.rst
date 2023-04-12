=================================
GWDALI Software
=================================

Installation
---------

Here is the console command:

.. code-block:: console

    $ pip install gwdali

Usage [example]
---------

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
	det0 = {'name':'CE','lon':45,'lat':45,'rot':0,'shape':90}
	# Einstein Telescope:
	det1 = {'name':'ET','lon':45,'lat':45,'rot':0,'shape':60}
	det2 = {'name':'ET','lon':45,'lat':45,'rot':120,'shape':60}
	det3 = {'name':'ET','lon':45,'lat':45,'rot':-120,'shape':60}

	#------------------------------------------------------
	# Setting Injections (Single detection)
	#------------------------------------------------------
	z = 0.1 # Redshift
	params = {}
	m1 = 1.5*(1+z)
	m2 = 1.5*(1+z)

	M = m1+m2 ; M2 = M*M 	# Total mass
	eta = m1*m2/M2 		# Symetric mass ratio
	Mc  = eta**(3./5)*M 	# Chirp Mass

	params['m1']  = m1
	params['m2']  = m2
	params['z']   = z
	params['RA']       = np.random.uniform(-180,180)
	params['Dec']      = (np.pi/2-np.arccos(np.random.uniform(-1,1)))*deg
	params['DL']       = cosmo.luminosity_distance(z).value/1.e3 # Gpc
	params['iota']     = np.random.uniform(0,np.pi) #np.pi/18 # 10 deg
	params['psi']      = np.random.uniform(-np.pi,np.pi)
	params['t_coal']   = 0
	params['phi_coal'] = 0
	params['sx1'] = 0 
	params['sy1'] = 0
	params['sz1'] = 0
	params['sx2'] = 0
	params['sy2'] = 0
	params['sz2'] = 0

	#----------------------------------------------------------------------
	# "approximant" options: 
	#		[Leading_Order, TaylorF2, TaylorF2_lal, IMRPhenomP, IMRPhenomD]
	#----------------------------------------------------------------------
	# "dali_method" options:
	#		[Fisher, Fisher_Sampling, Doublet, Triplet]
	#----------------------------------------------------------------------
	res = gw.GWDALI( Detection_Dict = params, 
			 FreeParams     = FreeParams, 
			 detectors      = [det0,det1,det2,det3], # Einstein Telescope + Cosmic Explorer
			 approximant    = "TaylorF2",
			 dali_method    = 'Fisher',
			 sampler_method = 'nestle', # Same as Bilby sampling method
			 save_fisher    = False,
			 save_cov       = False,
			 plot_corner    = False,
			 save_samples   = False,
			 hide_info      = True,
			 index          = 1,
			 r_cond		= 1.e-4,
			 npoints=300) # points for "nested sampling" or steps/walkers for "MCMC"

	Samples = res['Samples']
	Fisher  = res['Fisher']
	CovFish = res['CovFisher']
	Cov     = res['Covariance']
	Rec	= res['Recovery']
	Err     = res['Error']
	SNR     = res['SNR']

=================================  
API
=================================

.. py:function:: GWDALI.GWDALI(Detection_Dict, FreeParams, detectors, approximant, dali_method, samplet_method, save_fisher, save_cov, plot_corner, save_samples, hide_info, index, r_cond, npoints)

	Return GW samples, Fisher and covariance matrix, parameters uncertainties, parameters recovered and signal to noise ratio (SNR).

	:param Detection_Dict: A dictinary of GW parameters;
	:param FreeParams: list
	:param detectors: list
	:param approximant: string
	:param dali_method: string
	:param sampler_method: string
	:param save_fisher: boolean
	:param save_cov: boolean
	:param plot_corner: boolean
	:param save_samples: boolean
	:param hide_info: boolean
	:param index: integer
	:param r_cond: float
	:param npoints: integer
	:type Detection_Dict: dict

=================================  
About the Author
=================================

- Josiel Mendonça Soares de Souza (https://github.com/jmsdsouzaPhD)
- PhD in Physics by Universidade Federal do Rio Grande do Norte, Brazil
- Research Field: Gravitation, Cosmology and Gravitational Waves

=================================
License
=================================

MIT License

