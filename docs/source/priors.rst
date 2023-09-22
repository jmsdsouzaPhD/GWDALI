=================================  
Priors
=================================

.. code:: python

	Dict = {}
	DL0 = Detection_Dict['DL']
	d1 = cosmo.luminosity_distance(0.001).value / 1.e3 # Gpc
	d2 = cosmo.luminosity_distance(5.0).value / 1.e3   # Gpc
	#----------------------------#----------------------------#----------------------------
	Dict['DL']    = bilby.core.prior.Interped(name='DL',xx=xx,yy=yy,minimum=d1, maximum=d2)
	Dict['iota']  = bilby.core.prior.Sine(name='iota', minimum=0, maximum=np.pi)
	Dict['psi']   = bilby.core.prior.Uniform(name='psi',minimum=0, maximum=np.pi)
	#----------------------------#----------------------------#----------------------------
	Dict['alpha'] = bilby.core.prior.Uniform(name='alpha',minimum=-np.pi, maximum=np.pi)
	Dict['beta']  = bilby.core.prior.Sine(name='beta',minimum=0, maximum=np.pi)
	#----------------------------#----------------------------#----------------------------
	Dict['RA']    = bilby.core.prior.Uniform(name='RA',minimum=-180, maximum=180) # deg unit
	Dict['Dec']   = bilby.core.prior.Cosine(name='Dec',minimum=-np.pi/2, maximum=np.pi/2) # rad unit
	#----------------------------#----------------------------#----------------------------
	Dict['m1']    = bilby.core.prior.Uniform(name='m1',minimum=0.1, maximum=100)
	Dict['m2']    = bilby.core.prior.Uniform(name='m2',minimum=0.1, maximum=100)
	#----------------------------#----------------------------#----------------------------
	Dict['Mc']    = bilby.core.prior.Uniform(name='Mc',minimum=0.1, maximum=100)
	Dict['eta']   = bilby.core.prior.Uniform(name='eta',minimum=1.e-3, maximum=1./4)
	Dict['q']     = bilby.core.prior.Uniform(name='q',minimum=1.e-3, maximum=1.0)
	#----------------------------#----------------------------#----------------------------
	Dict['sx1']    = bilby.core.prior.Uniform(name='sx1',minimum=0, maximum=1.0)
	Dict['sy1']    = bilby.core.prior.Uniform(name='sy1',minimum=0, maximum=1.0)
	Dict['sz1']    = bilby.core.prior.Uniform(name='sz1',minimum=0, maximum=1.0)
	Dict['sx2']    = bilby.core.prior.Uniform(name='sx2',minimum=0, maximum=1.0)
	Dict['sy2']    = bilby.core.prior.Uniform(name='sy2',minimum=0, maximum=1.0)
	Dict['sz2']    = bilby.core.prior.Uniform(name='sz2',minimum=0, maximum=1.0)
	#----------------------------#----------------------------#----------------------------
	Dict['phi_coal']  = bilby.core.prior.Uniform(name='phi_coal',minimum=0, maximum=2*np.pi)
	Dict['t_coal']    = bilby.core.prior.Uniform(name='t_coal',minimum=0, maximum=3600) # 1 hour
	#----------------------------#----------------------------#----------------------------