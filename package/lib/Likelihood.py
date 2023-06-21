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
			if(fp=='Dec'): params[fp] *= 180/np.pi
			
		loglike = 0 ; ndet = 0
		detectors, approximant = self.DetAp 
		for det in detectors:
			h = gwfunc.Signal(params, det, approximant)
			loglike += -gwfunc.ScalarProduct(det['freq'],det['Sn'],self.Data[ndet]-h,self.Data[ndet]-h)
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

		Fisher_vec = Fisher.ravel() # convert 2D matrix to 1D array
		Doublet3, Doublet4 = Doublet
		Triplet4, Triplet5, Triplet6 = Triplet

		# dX_ijk:  Tensor rank 3
		# dX_ijkl: Tensor rank 4

		dX2 = np.zeros(Np*Np)
		if(self.dali_method in ['Doublet','Triplet']):
			dX3 = np.zeros(Np**3)
			dX4 = np.zeros(Np**4)
			if(self.dali_method == 'Triplet'):
				dX5 = np.zeros(Np**5)
				dX6 = np.zeros(Np**6)

		#-------------------------------------------------------------
		# For Fisher
		#-------------------------------------------------------------
		for i in range(Np):
			for j in range(Np):
				idx2 = Np*i + j
				dX2[idx2] = dT[i]*dT[j]
				#-------------------------------------------------------------
				# For Doublet
				#-------------------------------------------------------------
				if(self.dali_method in ['Doublet','Triplet']):
					for k in range(Np):
						idx3 = i*Np**2 + j*Np + k
						dX3[idx3] = dT[i]*dT[j]*dT[k]
						for l in range(Np):
							idx4 = i*Np**3 + j*Np**2 + k*Np + l
							dX4[idx4] = dT[i]*dT[j]*dT[k]*dT[l]
							#-------------------------------------------------------------
							# For Triplet
							#-------------------------------------------------------------
							if(self.dali_method == 'Triplet'):
								for m in range(Np):
									idx5 = i*Np**4 + j*Np**3 + k*Np**2 + l*Np + m
									dX5[idx5] = dT[i]*dT[j]*dT[k]*dT[l]*dT[m]
									for n in range(Np):
										idx6 = i*Np**5 + j*Np**4 + k*Np**3 + l*Np**2 + m*Np + n
										dX6[idx6] = dT[i]*dT[j]*dT[k]*dT[l]*dT[m]*dT[n]
		#-------------------------------------------------------------

		loglike = - sum(Fisher_vec*dX2)/2.
		if(self.dali_method in ['Doublet','Triplet']):
			loglike -= (  sum(Doublet3*dX3)/2. + sum(Doublet4*dX4)/8. )
			if(self.dali_method == 'Triplet'):
				loglike -= ( sum(Triplet4*dX4)/6. + sum(Triplet5*dX5)/12. + sum(Triplet6*dX6)/72. )

		if(np.isnan(loglike) or np.isnan(loglike)): loglike = -1.e10
		return loglike

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

def get_posterior(FreeParams, Theta0, Detection_Dict, GwData, approximant, Detectors, Tensors, dali_method, sampler_method, npoints,new_priors):
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