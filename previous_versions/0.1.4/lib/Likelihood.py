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
		self.h_data       = args[1]
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

		det_a = detectors[0]
		for det_b in detectors:
			dets = [det_a, det_b]
			h_model = gwfunc.Signal(params, dets, approximant)#*ExpT

			diff = self.h_data[ndet]-h_model
			loglike -= 0.5*gwfunc.ScalarProduct(det_b['freq'],det_b['Sn'],diff,diff)
			ndet += 1
			
		if(np.isnan(loglike)): loglike = -np.inf
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
		Db12, Db22 = Doublet
		Tp13, Tp23, Tp33 = Triplet

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
		# logL = logL0 - (1/2)Fisher * dTheta_ij - [ (1/2)Db12 * dTheta_ijk  + (1/8) * dTheta_ijkl ]
		#        - [ (1/6)Tp13 * dTheta_ijkl + (1/12) * dTheta_ijklm + (1/72) * dTheta_ijklmn ] + ...
		
		loglike = -np.sum(Fij*dX2)/2
		loglike -= np.sum(Db12*dX3)/2 + np.sum(Db22*dX4)/8
		loglike -= np.sum(Tp13*dX4)/6 + np.sum(Tp23*dX5)/12 + np.sum(Tp33*dX6)/72

		if(np.isnan(loglike)): loglike = -np.inf
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
# Prior(1/DL): dVc/dz * dL^2
# p(x)dx = p(x)|dx/dy| dy -> p(y) = p(x)|dx/dy|
#---------------------------------------- 
XidL = 1./XdL ; idxs = np.argsort(XidL)
XidL = XidL[idxs]
YidL = YdL[idxs] * XdL[idxs]**2
YidL = YidL/sum(YidL)
#----------------------------------------
# Prior(logDL): dVc/dz * dL
#---------------------------------------- 
XlogDL = np.log(XdL) ; idxs = np.argsort(XlogDL)
XlogDL = XlogDL[idxs]
YlogDL = YdL[idxs] * XdL[idxs]
YlogDL = YlogDL/sum(YlogDL)
#----------------------------------------
d1, d2 = XdL[0] , XdL[-1]
id1, id2 = XidL[0] , XidL[-1]
lD1, lD2 = XlogDL[0] , XlogDL[-1]
#----------------------------------------
# Pior(Dec): cos(Dec[degrees])
Xdec = np.linspace(-90,90,1000)
Ydec = np.cos(Xdec*np.pi/180)
Ydec/=sum(Ydec)
#----------------------------------------

Priors_std = {}
Priors_std['DL']	   = bilby.core.prior.Interped(name='DL',xx=XdL,yy=YdL,minimum=d1, maximum=d2) # Gpc
Priors_std['inv_dL']   = bilby.core.prior.Interped(name='inv_dL',xx=XidL,yy=YidL,minimum=id1, maximum=id2) # 1./Gpc
Priors_std['ln_dL']    = bilby.core.prior.Interped(name='ln_dL',xx=XlogDL,yy=YlogDL,minimum=lD1, maximum=lD2) # ln(Gpc)
Priors_std['iota']	   = bilby.core.prior.Sine(name='iota', minimum=0, maximum=np.pi) # radians
Priors_std['cos_iota'] = bilby.core.prior.Uniform(name='cos_iota',minimum=-1., maximum=1.) # radians
Priors_std['psi']      = bilby.core.prior.Uniform(name='psi',minimum=0, maximum=np.pi) # radians
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
#----------------------------
Priors_std['S1']		= bilby.core.prior.Uniform(name='S1',minimum=0, maximum=.9)
Priors_std['theta_1']	= bilby.core.prior.Sine(name='theta_1',minimum=0, maximum=np.pi)
Priors_std['phi_1']		= bilby.core.prior.Uniform(name='phi_1',minimum=0, maximum=2*np.pi)
Priors_std['S2']		= bilby.core.prior.Uniform(name='S2',minimum=0, maximum=.9)
Priors_std['theta_2']	= bilby.core.prior.Sine(name='theta_2',minimum=0, maximum=np.pi)
Priors_std['phi_2']		= bilby.core.prior.Uniform(name='phi_2',minimum=0, maximum=2*np.pi)
#----------------------------#----------------------------#----------------------------
Priors_std['phi_coal']  = bilby.core.prior.Uniform(name='phi_coal',minimum=0, maximum=2*np.pi) # radians
Priors_std['t_coal']    = bilby.core.prior.Uniform(name='t_coal',minimum=0, maximum=3600) # seconds
#----------------------------#----------------------------#----------------------------

def get_posterior(FreeParams, Theta0, Detection_Dict, GwData, approximant, Detectors, Tensors, dali_method, sampler_method, npoints,new_priors):
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

	outdir = 'outdir_%s_%d/' % (dali_method, int(np.random.uniform(0,1000))  )
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
