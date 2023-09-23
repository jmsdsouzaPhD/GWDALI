=================================  
Priors
=================================

.. code:: python

	import bilby
	import numpy as np
	from astropy.cosmology import FlatLambdaCDM
	cosmo = FlatLambdaCDM(70,0.3)

	Priors = {}

	#----------------------------------------
	# Prior(DL): dVc/dz
	#----------------------------------------
	x = np.linspace(1.e-3,10,1000)
	d = cosmo.luminosity_distance(x).value
	H = cosmo.H(x).value
	Rc = d/(1+x)
	priorD = 4*np.pi*Rc**2/H
	xx = d/1.e3
	yy = priorD/sum(priorD)
	#----------------------------#----------------------------#----------------------------
	d1 = cosmo.luminosity_distance(0.001).value / 1.e3 # Gpc
	d2 = cosmo.luminosity_distance(5.0).value / 1.e3   # Gpc
	Priors['DL']    = bilby.core.prior.Interped(name='DL',xx=xx,yy=yy,minimum=d1, maximum=d2)
	Priors['iota']  = bilby.core.prior.Sine(name='iota', minimum=0, maximum=np.pi)
	Priors['psi']   = bilby.core.prior.Uniform(name='psi',minimum=0, maximum=np.pi)
	#----------------------------#----------------------------#----------------------------
	Priors['alpha'] = bilby.core.prior.Uniform(name='alpha',minimum=-np.pi, maximum=np.pi)
	Priors['beta']  = bilby.core.prior.Sine(name='beta',minimum=0, maximum=np.pi)
	#----------------------------#----------------------------#----------------------------
	xx = np.linspace(-90,90,10000)
	yy = np.cos(xx*np.pi/180)
	yy/=sum(yy)
	Priors['RA']    = bilby.core.prior.Uniform(name='RA',minimum=-180, maximum=180) # deg unit
	Priors['Dec']   = bilby.core.prior.Interped(name='Dec',xx=xx,yy=yy,minimum=-90, maximum=90) # deg unit
	#----------------------------#----------------------------#----------------------------
	Priors['m1']    = bilby.core.prior.Uniform(name='m1',minimum=0.1, maximum=100)
	Priors['m2']    = bilby.core.prior.Uniform(name='m2',minimum=0.1, maximum=100)
	#----------------------------#----------------------------#----------------------------
	Priors['Mc']    = bilby.core.prior.Uniform(name='Mc',minimum=0.1, maximum=100)
	Priors['eta']   = bilby.core.prior.Uniform(name='eta',minimum=1.e-3, maximum=1./4)
	Priors['q']     = bilby.core.prior.Uniform(name='q',minimum=1.e-3, maximum=1.0)
	#----------------------------#----------------------------#----------------------------
	Priors['sx1']    = bilby.core.prior.Uniform(name='sx1',minimum=0, maximum=1.0)
	Priors['sy1']    = bilby.core.prior.Uniform(name='sy1',minimum=0, maximum=1.0)
	Priors['sz1']    = bilby.core.prior.Uniform(name='sz1',minimum=0, maximum=1.0)
	Priors['sx2']    = bilby.core.prior.Uniform(name='sx2',minimum=0, maximum=1.0)
	Priors['sy2']    = bilby.core.prior.Uniform(name='sy2',minimum=0, maximum=1.0)
	Priors['sz2']    = bilby.core.prior.Uniform(name='sz2',minimum=0, maximum=1.0)
	#----------------------------#----------------------------#----------------------------
	Priors['phi_coal']  = bilby.core.prior.Uniform(name='phi_coal',minimum=0, maximum=2*np.pi)
	Priors['t_coal']    = bilby.core.prior.Uniform(name='t_coal',minimum=0, maximum=3600) # 1 hour
	#----------------------------#----------------------------#----------------------------

.. figure:: ./fig_priors.png
   :alt: priors
   :align: center