from __future__ import division
import numpy as np
import bilby, os, shutil, sys
from time import time as now

import GWDALI.lib.Diff_Signal_Tensors as gwfunc
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
		for f in args[0]: aux[f] = None
		super().__init__(aux)

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
			
		if(np.isnan(loglike) or np.isinf(loglike)): loglike = -1.e9
		return loglike
#--------------------------------------------------

class DALI_likelihood(bilby.Likelihood):
	def __init__(self,args):
		self.FreeParams  = args[0]
		self.Theta0      = args[1]
		self.Tensors     = args[2]
		self.dali_method = args[3]
		aux = {}
		for f in args[0]: aux[f] = None
		super().__init__(aux)

	def log_likelihood(self):
		dX = [] ; i = 0
		Np = len(self.FreeParams)
		for fp in self.FreeParams:
			Theta  = self.parameters[fp]
			Theta0 = self.Theta0[i] ; i += 1
			dX.append(Theta-Theta0) # Fixed!!!

		Fij , Doublet , Triplet = self.Tensors
		Db3, Db4 = Doublet
		Tp4, Tp5, Tp6 = Triplet

		dX  = np.array(dX)
		dX2 = np.outer(dX,dX)
		dX3 = np.outer(dX2,dX)
		dX4 = np.outer(dX3,dX)
		dX5 = np.outer(dX4,dX)
		dX6 = np.outer(dX5,dX)

		dX2 = np.ravel(dX2)
		dX3 = np.ravel(dX3)
		dX4 = np.ravel(dX4)
		dX5 = np.ravel(dX5)
		dX6 = np.ravel(dX6)

		# From arXiv:2203.02670
		# logL = logL0 - (1/2)Fisher * dTheta_ij - [ (1/2)Db3 * dTheta_ijk  + (1/8) * dTheta_ijkl ]
		#        - [ (1/6)Tp4 * dTheta_ijkl + (1/12) * dTheta_ijklm + (1/72) * dTheta_ijklmn ] + ...
		
		loglike = -np.sum(Fij*dX2)/2
		loglike -= np.sum(Db3*dX3)/2 + np.sum(Db4*dX4)/8
		loglike -= np.sum(Tp4*dX4)/6 + np.sum(Tp5*dX5)/12 + np.sum(Tp6*dX6)/72

		if(np.isnan(loglike) or np.isnan(loglike)): return -1.e9
		else: return loglike

#----------------------------------------
# Prior(DL): dVc/dz
#---------------------------------------- 
z      = np.linspace(1.e-3,10,1000) ; dz = z[1]-z[0]# redshift
Vc     = cosmo.comoving_volume(z).value
priorD = np.diff(Vc)/np.diff(z)
z   = z[1:]-dz/2
XdL = cosmo.luminosity_distance(z).value / 1.e3 # dL in Gpc
YdL = priorD/sum(priorD)
#----------------------------------------
d1 = cosmo.luminosity_distance(1.e-3).value / 1.e3 # Gpc
d2 = cosmo.luminosity_distance(5.0).value / 1.e3   # Gpc
#----------------------------------------
# Pior(Dec): cos(Dec[degrees])
Xdec = np.linspace(-90,90,1000)
Ydec = np.cos(Xdec*np.pi/180)
Ydec/=sum(Ydec)
#----------------------------------------

def get_posterior(FreeParams, Theta0, Detection_Dict, GwData, approximant, Detectors, Tensors, dali_method, sampler_method, npoints,new_priors):
	Priors_std = {}
	Priors_std['DL']   = bilby.core.prior.Interped(name='DL',xx=XdL,yy=YdL,minimum=d1, maximum=d2) # Gpc
	Priors_std['iota']  = bilby.core.prior.Sine(name='iota', minimum=0, maximum=np.pi) # radians
	Priors_std['psi']   = bilby.core.prior.Uniform(name='psi',minimum=0, maximum=np.pi) # radians
	#----------------------------#----------------------------#----------------------------
	Priors_std['alpha'] = bilby.core.prior.Uniform(name='alpha',minimum=-np.pi, maximum=np.pi) # radians
	Priors_std['beta']  = bilby.core.prior.Sine(name='beta',minimum=0, maximum=np.pi) # radians
	#----------------------------#----------------------------#----------------------------
	Priors_std['RA']    = bilby.core.prior.Uniform(name='RA',minimum=-180, maximum=180) # degrees
	Priors_std['Dec']   = bilby.core.prior.Interped(name='Dec',xx=Xdec,yy=Ydec,minimum=-90, maximum=90) # degrees
	#----------------------------#----------------------------#----------------------------
	Priors_std['m1']    = bilby.core.prior.Uniform(name='m1',minimum=0.1, maximum=100) # solar mass
	Priors_std['m2']    = bilby.core.prior.Uniform(name='m2',minimum=0.1, maximum=100) # solar mass
	#----------------------------#----------------------------#----------------------------
	Priors_std['Mc']    = bilby.core.prior.Uniform(name='Mc',minimum=0.1, maximum=100) # solar mass
	Priors_std['eta']   = bilby.core.prior.Uniform(name='eta',minimum=1.e-3, maximum=1./4)
	Priors_std['q']     = bilby.core.prior.Uniform(name='q',minimum=1.e-3, maximum=1.0)
	#----------------------------#----------------------------#----------------------------
	Priors_std['sx1']    = bilby.core.prior.Uniform(name='sx1',minimum=0, maximum=1.0)
	Priors_std['sy1']    = bilby.core.prior.Uniform(name='sy1',minimum=0, maximum=1.0)
	Priors_std['sz1']    = bilby.core.prior.Uniform(name='sz1',minimum=0, maximum=1.0)
	Priors_std['sx2']    = bilby.core.prior.Uniform(name='sx2',minimum=0, maximum=1.0)
	Priors_std['sy2']    = bilby.core.prior.Uniform(name='sy2',minimum=0, maximum=1.0)
	Priors_std['sz2']    = bilby.core.prior.Uniform(name='sz2',minimum=0, maximum=1.0)
	#----------------------------#----------------------------#----------------------------
	Priors_std['phi_coal']  = bilby.core.prior.Uniform(name='phi_coal',minimum=0, maximum=2*np.pi) # radians
	Priors_std['t_coal']    = bilby.core.prior.Uniform(name='t_coal',minimum=0, maximum=3600) # seconds
	#----------------------------#----------------------------#----------------------------

	Priors = {}
	for fp in FreeParams: Priors[fp] = Priors_std[fp]
	if(new_priors!=None):
		for fp in new_priors.keys():
			try:
				Xval, Yval = new_priors[fp]
				Yval /= sum(Yval)
				Priors[fp] = bilby.core.prior.Interped(name=fp,xx=Xval,yy=Yval,minimum=min(Xval),maximum=max(Xval))
			except:
				print('\n>> Invalid key/argment in new_priors dictionay!\n')
				quit()

	if(dali_method=='Standard'):
		likelihood = GW_likelihood([FreeParams, GwData, Detection_Dict, [Detectors,approximant] ])
	else:
		Fisher, Doublet, Triplet = Tensors
		Fisher = np.ravel(Fisher)
		Doublet = [np.ravel(Doublet[i]) for i in range(2)]
		Triplet = [np.ravel(Triplet[i]) for i in range(3)]
		Tensors = [Fisher,Doublet,Triplet]
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
							nburn=int(0.2*npoints), outdir=outdir)
	shutil.rmtree(outdir)

	#----------------------------
	M = [res.posterior[fp] for fp in FreeParams]
	return np.transpose(np.array(M))
