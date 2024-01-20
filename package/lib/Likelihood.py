from __future__ import division
import numpy as np
import bilby, os, shutil, sys
from time import time as now

import GWDALI.lib.Derivatives_Tensors as gwfunc
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(70,0.3)

#sys.stdout = open(os.devnull, 'w')
bilby.core.utils.setup_logger(log_level=0)	# hide init bilby output in screen
bilby.core.utils.setup_logger(log_level=10)
bilby.core.utils.setup_logger(log_level=20)
bilby.core.utils.setup_logger(log_level=30)
bilby.core.utils.setup_logger(log_level=40)
bilby.core.utils.setup_logger(log_level=50)

#--------------------------------------------------
class GW_likelihood(bilby.Likelihood):
	def __init__(self,args):
		aux = {}
		self.FreeParams = args[0]
		self.Data       = args[1]
		self.GWparams   = args[2]
		self.DetAp		= args[3]
		for f in args[0]: aux[f] = 0
		bilby.Likelihood.__init__(self, aux)
		pass

	def log_likelihood(self):
		Np = len(self.FreeParams)
		params = self.GWparams.copy()
		for fp in self.FreeParams:
			params[fp] = self.parameters[fp]
			
		loglike = 0 ; ndet = 0
		detectors, approximant = self.DetAp 

		for det in detectors:
			h = gwfunc.Signal(params, det, approximant)
			diff = self.Data[ndet]-h
			loglike -= 0.5*gwfunc.ScalarProduct(det['freq'],det['Sn'],diff,diff)
			ndet += 1
			
		if(np.isnan(loglike) or np.isinf(loglike)): loglike = -1.e10
		return loglike
#--------------------------------------------------

class DALI_likelihood(bilby.Likelihood):
	def __init__(self,args):
		self.FreeParams  = args[0]
		self.Theta0      = args[1]
		self.Tensors     = args[2]
		self.dali_method = args[3]

		aux = {}
		for f in args[0]:
			aux[f] = 0

		self.it = 0
		bilby.Likelihood.__init__(self, aux)
		pass

	def log_likelihood(self):
		Theta = [] ; Np = len(self.FreeParams)
		for fp in self.FreeParams:
			prm = self.parameters[fp]
			if(fp=='Dec'): prm *= 180/np.pi
			Theta.append(prm)

		Theta  = np.array(Theta)
		dT = np.array( self.Theta0 - Theta )

		Fisher , Doublet , Triplet = self.Tensors

		Fisher_vec = Fisher.ravel() # ravel() used to convert n-dimensional matrix to 1D array
		Doublet3, Doublet4 = Doublet
		Triplet4, Triplet5, Triplet6 = Triplet

		dX2 = np.outer(dT,dT)
		dX3 = np.outer(dX2,dT)
		dX4 = np.outer(dX3,dT)
		dX5 = np.outer(dX4,dT)
		dX6 = np.outer(dX5,dT)

		dX2 = dX2.ravel()
		dX3 = dX3.ravel()
		dX4 = dX4.ravel()
		dX5 = dX5.ravel()
		dX6 = dX6.ravel()

		# From arXiv:2203.02670
		# logL = logL0 - (1/2)Fisher * dTheta_ij - [ (1/2)Doublet3 * dTheta_ijk  + (1/8) * dTheta_ijkl ]
		#        - [ (1/6)Triplet4 * dTheta_ijkl + (1/12) * dTheta_ijklm + (1/72) * dTheta_ijklmn ] + ...

		loglike = - sum(Fisher_vec*dX2)/2.
		if(self.dali_method in ['Doublet','Triplet']):
			loglike -= (  sum(Doublet3*dX3)/2. + sum(Doublet4*dX4)/8. )
			if(self.dali_method == 'Triplet'):
				loglike -= ( sum(Triplet4*dX4)/6. + sum(Triplet5*dX5)/12. + sum(Triplet6*dX6)/72. )

		if(np.isnan(loglike) or np.isnan(loglike)): loglike = -1.e10
		return loglike

def get_posterior(FreeParams, Theta0, Detection_Dict, GwData, approximant, Detectors, Tensors, dali_method, sampler_method, npoints,new_priors):
	Dict = {}
	d1 = cosmo.luminosity_distance(0.001).value / 1.e3 # Gpc
	d2 = cosmo.luminosity_distance(5.0).value / 1.e3   # Gpc
	#----------------------------#----------------------------#----------------------------
	
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
	#----------------------------------------

	Dict['DL']    = bilby.core.prior.Interped(name='DL',xx=xx,yy=yy,minimum=d1, maximum=d2)
	Dict['iota']  = bilby.core.prior.Sine(name='iota', minimum=0, maximum=np.pi)
	Dict['psi']   = bilby.core.prior.Uniform(name='psi',minimum=0, maximum=np.pi)
	#----------------------------#----------------------------#----------------------------
	Dict['alpha'] = bilby.core.prior.Uniform(name='alpha',minimum=-np.pi, maximum=np.pi)
	Dict['beta']  = bilby.core.prior.Sine(name='beta',minimum=0, maximum=np.pi)
	#----------------------------#----------------------------#----------------------------
	xx = np.linspace(-90,90,1000)
	yy = np.cos(xx*np.pi/180)
	yy/=sum(yy)
	Dict['RA']    = bilby.core.prior.Uniform(name='RA',minimum=-180, maximum=180) # deg unit
	Dict['Dec']   = bilby.core.prior.Interped(name='Dec',xx=xx,yy=yy,minimum=-90, maximum=90) # deg unit
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

	Priors = {}
	for fp in FreeParams:
		if(new_priors == None):
			Priors[fp] = Dict[fp]
		elif(not (fp in new_priors.keys()) ):
			Priors[fp] = Dict[fp]
		else:
			try:
				x, y = new_priors[fp]
				y /= sum(y)
				Priors[fp] = bilby.core.prior.Interped(name=fp,xx=x,yy=y,minimum=min(x),maximum=max(x))
			except:
				print('\n>> Invalid key/argment in new_priors dictionay!\n')
				quit()

	if(dali_method=='Standard'):
		likelihood = GW_likelihood([FreeParams, GwData, Detection_Dict, [Detectors,approximant] ])
	else:
		likelihood = DALI_likelihood([FreeParams, Theta0, Tensors, dali_method])

	outdir = 'outdir_%s_%d/' % (dali_method, int(np.random.uniform(0,300))  )
	if(os.path.isdir(outdir)):
		shutil.rmtree(outdir)

	#----------------------------
	# Running Sampler
	#----------------------------
	print("\t >> Method:", dali_method)
	res = bilby.run_sampler(likelihood=likelihood, priors=Priors,
							sampler=sampler_method, npoints=npoints, 
							nsteps=npoints, nwalkers=npoints, 
							nburn=int(0.8*npoints), outdir=outdir)
	shutil.rmtree(outdir)

	#----------------------------

	M = []
	for fp in FreeParams:
		M.append(res.posterior[fp])
	M = np.transpose(np.array(M))

	return M