import numpy as np
import bilby, os
import jax.numpy as jnp
from jax import jit, vmap

from . import Waveforms as wf
from . import Tensors as Tns
from . import GridLike as grd

import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from scipy.ndimage import gaussian_filter as gf
from scipy.interpolate import interp1d
from tqdm import trange
from time import time as now
import traceback

bilby.core.utils.setup_logger(log_level=0)	# hide init bilby output in screen
bilby.core.utils.setup_logger(log_level=10)
bilby.core.utils.setup_logger(log_level=20)
bilby.core.utils.setup_logger(log_level=30)
bilby.core.utils.setup_logger(log_level=40)
bilby.core.utils.setup_logger(log_level=50)

Kappa = 7.776922780044755e-23

R_earth_equatiorial = 6.3781e6
R_earth_polar = 6.3568e6
R_earth = 6.378e6 # meters
c = 299792458. # speed of light

#--------------------------------------------------------------------
# -------------------------- GW Parameters --------------------------
#--------------------------------------------------------------------

try:
	from jax.scipy.integrate import trapezoid
except:
	print("Fail on loading jax.scypi.integrate.trapezoid")
	@jit
	def trapezoid(y,x):
		dx = jnp.diff(x)
		return jnp.sum( 0.5 * (y[1:]+y[:-1])*dx )

class Exact_likelihood(bilby.Likelihood): 
	def __init__(self,args):
		self.Data = args[0]
		self.Dets = args[1]
		self.Prms = args[2]
		self.GwSignal = args[3]
		self.FreeParams = args[4]
		self.vec_time = []
		aux = {}
		for x in self.FreeParams: aux[x] = None
		self.cont = 0
		super().__init__(aux)

	def log_likelihood(self):
		try:
			time_init = now()
			p = self.Prms.copy()
			for x in self.FreeParams:
				p[x] = self.parameters[x]

			gw_values = list(p.values())
			loglike = 0
			det_conf_ref = [self.Dets[0][x] for x in "lon,lat,rot,shape".split(',')]
			for n, det in enumerate(self.Dets):
				H = self.Data[n]
				Sn, freq = Tns.get_Sn(det['name'])
				det_conf = [ [det[x] for x in "lon,lat,rot,shape".split(',')] , det_conf_ref ]

				h  = self.GwSignal(gw_values, det_conf, freq)

				dH = H-h
				loglike -= 0.5*Tns.ScalarProd(dH,dH,Sn,freq)
			self.vec_time.append( now() - time_init )
			self.cont += 1
			if(np.isnan(loglike)): return -1.e10

		except Exception as e:
			loglike = -1.e10 ; self.cont += 1
			print(f"[{self.cont}] Error in likelihood evaluation!\n{e}\n")
			traceback.print_exc()
			print()

			print("likelihood Fail!") ; quit()
		return loglike

	def get_times(self):
		return self.vec_time

@jit
def compute_logL1(Fisher,dT_vec):
	dT,dT2,dT3,dT4,dT5,dT6 = dT_vec
	return -0.5*jnp.sum( Fisher*dT2 )

@jit
def compute_logL2(Fisher,Db,dT_vec):
	dT,dT2,dT3,dT4,dT5,dT6 = dT_vec
	Db12, Db22 = Db
	logL1 = -0.5*jnp.sum( Fisher*dT2 )
	logL2 = -0.5*jnp.sum( Db12*dT3 ) - 0.125*jnp.sum( Db22*dT4 )
	return logL1 + logL2

@jit
def compute_logL3(Fisher,Db,Tp,dT_vec):
	dT,dT2,dT3,dT4,dT5,dT6 = dT_vec
	Db12, Db22 = Db
	Tp13, Tp23, Tp33 = Tp

	logL1 = -0.5*jnp.sum( Fisher*dT2 )
	logL2 = -0.5*jnp.sum( Db12*dT3 ) - 0.125*jnp.sum( Db22*dT4 )
	logL3 = -(1./6)*jnp.sum( Tp13*dT4 ) - (1./12)*jnp.sum( Tp23*dT5 ) - (1./72)*jnp.sum( Tp33*dT6 )
	return logL1 + logL2 + logL3

class DALI_likelihood(bilby.Likelihood):
	# Fihser/DALI Likelihoods
	def __init__(self,args):
		self.gwprms     = args[0]
		self.Tensors    = args[1]
		self.FreeParams = args[2]
		self.Method     = args[3]
		self.vec_time   = []
		aux = {}
		for fp in self.FreeParams: aux[fp] = None
		super().__init__(aux)

	def log_likelihood(self):
		time_init = now()

		dT  = jnp.array([self.parameters[fp] - self.gwprms[fp] for fp in self.FreeParams])
		dT2 = jnp.tensordot(dT,dT,axes=0)
		dT3 = jnp.tensordot(dT2,dT,axes=0)
		dT4 = jnp.tensordot(dT3,dT,axes=0)
		dT5 = jnp.tensordot(dT4,dT,axes=0)
		dT6 = jnp.tensordot(dT5,dT,axes=0)

		dT2 = np.ravel(dT2)
		dT3 = np.ravel(dT3)
		dT4 = np.ravel(dT4)
		dT5 = np.ravel(dT5)
		dT6 = np.ravel(dT6)

		Fisher, Doublet, Triplet = self.Tensors
		Db12, Db22       = Doublet
		Tp13, Tp23, Tp33 = Triplet

		#======================================
		# Speedup likelihood with @jit
		#======================================
		dT_vec = dT,dT2,dT3,dT4,dT5,dT6
		if self.Method == "Fisher": 
			loglike = compute_logL1(Fisher,dT_vec)
		elif self.Method == "Doublet": 
			loglike = compute_logL2(Fisher,Doublet,dT_vec)
		elif self.Method == "Triplet": 
			loglike = compute_logL3(Fisher,Doublet,Triplet,dT_vec)
		else:
			print("Method Invalid")
			quit()

		self.vec_time.append( now() - time_init )
		if(np.isnan(loglike)): return -1.e10
		else: return loglike

	def get_times(self):
		return self.vec_time

# Prior to Dec:
xx = np.linspace(-90,90,1000)
yy = np.cos(xx*np.pi/180)
yy/=np.sum(yy)

def Normed(x,y):
	ii = np.argsort(x)
	x, y = x[ii], y[ii]
	return x, y/trapezoid(y,x)

#=======================================
# Prior to dL:
#=======================================
DL_TRANSFORMS = {"dL":         dict( Y = lambda x: x , X = lambda y: y, J = lambda x: 1.0 ),
				 "lnDL":      dict( Y = lambda x: np.log(x) , X = lambda y: np.exp(y), J = lambda x: x ),
				 "inv_dL":     dict( Y = lambda x: 1./x , X = lambda y: 1./y , J = lambda x: x**2),
				 "inv_dL2":    dict( Y = lambda x: 1./x**2 , X = lambda y: 1./y**.5 , J = lambda x: .5*x**3),
				 "inv_sqrtdL": dict( Y = lambda x: 1./x**.5 , X = lambda y: 1./y**2 , J = lambda x: 2.*x**1.5),
				 "inv_lnDL":   dict( Y = lambda x: 1./np.log(x) , X = lambda y: np.exp(1./y) , J = lambda x: x*np.log(x)**2),
				}

H0, Om = 70,0.3
cKm = c/1.e3
cosmo = FlatLambdaCDM(H0,Om)	

def P_dL(y,key):

	K = (H0/cKm)*1.e3
	T = DL_TRANSFORMS[key]
	X = T["X"](y)

	z = K*X
	DL = cosmo.luminosity_distance(z).value/1.e3

	func_z = interp1d(DL,z,fill_value=0,bounds_error=False)

	Vc = cosmo.comoving_volume(z).value
	P_dL = np.diff(Vc) / np.diff(DL)
	DL = .5*(DL[1:]+DL[:-1])
	
	zi = func_z(DL)
	sfr = (1.+zi)**(2.7) / ( 1. + ((1.+zi)/2.9)**5.6)
	P_dL *= sfr/(1+zi)

	u = T["Y"](DL)
	v = P_dL * T["J"](DL)

	idxs = np.argsort(u)
	u = u[idxs]
	v = v[idxs]
	
	return u, v/trapezoid(v,u)

N = 10000
y1 = np.linspace(1.e-3,100,N)

ys = { "dL": y1,
	   "inv_dL": 1./y1,
	   "inv_dL2": 1./y1**2,
	   "inv_sqrtdL": 1./y1**.5,
	   "inv_lnDL": 1./np.log(y1),
	   "lnDL": np.log(y1)
	  }

PriorsX = {}

for key, y in ys.items():
	u, v = P_dL(y, key)
	PriorsX[key] = bilby.core.prior.Interped(name=key,xx=u,yy=v,minimum=u.min(),maximum=u.max())

#=======================================
eps = 1.e-2
q = np.logspace(-2,np.log10(1-eps),10000)
eta = q/(1+q)**2
inv_eta = 1./eta
ln_eta = np.log(eta)
deltaM = np.sqrt(1-4*eta)

Pq = q**2 ; Pq/=np.sum(Pq)
P_eta = Pq * (1.+q)**3/(1.-q) ; P_eta/=np.sum(P_eta)
P_invEta = P_eta * eta**2 ; P_invEta/=np.sum(P_invEta)
P_lnEta = P_eta * eta ; P_lnEta/=np.sum(P_lnEta)
P_deltaM = P_eta * .5*deltaM

Mc = np.linspace(1,100,10000)
P_Mc = np.ones_like(Mc)
ln_Mc = np.log(Mc)
P = Mc*P_Mc
P_ln_Mc = P/np.sum(P)

PriorsX['m1'] = bilby.core.prior.Uniform(name='m1',minimum=1,maximum=100)
PriorsX['m2'] = bilby.core.prior.Uniform(name='m2',minimum=1,maximum=100)
PriorsX['M']  = bilby.core.prior.Uniform(name='M',minimum=1,maximum=200)

PriorsX['Mc']    = bilby.core.prior.Uniform(name='Mc',minimum=1,maximum=100)
PriorsX['ln_Mc'] = bilby.core.prior.Interped(name='ln_Mc',xx=ln_Mc,yy=P_ln_Mc,minimum=ln_Mc.min(),maximum=ln_Mc.max())

PriorsX['q']       = bilby.core.prior.Interped(name='q',xx=q,yy=Pq, minimum=q.min(),maximum=q.max())
PriorsX['eta']     = bilby.core.prior.Interped(name='eta',xx=eta,yy=P_eta, minimum=eta.min(),maximum=eta.max())
PriorsX['inv_eta'] = bilby.core.prior.Interped(name='inv_eta',xx=inv_eta,yy=P_invEta,minimum=inv_eta.min(),maximum=inv_eta.max())
PriorsX['ln_eta']  = bilby.core.prior.Interped(name='ln_eta',xx=ln_eta,yy=P_lnEta,minimum=ln_eta.min(),maximum=ln_eta.max())
PriorsX['deltaM']  = bilby.core.prior.Interped(name='deltaM',xx=deltaM,yy=P_deltaM,minimum=deltaM.min(),maximum=deltaM.max())

PriorsX['RA']       = bilby.core.prior.Uniform(name='RA',minimum=0,maximum=360)
PriorsX['Dec']      = bilby.core.prior.Interped(name='Dec',xx=xx,yy=yy,minimum=-90,maximum=90)
PriorsX['iota']     = bilby.core.prior.Sine(name='iota',minimum=0,maximum=np.pi)
PriorsX['cos_iota'] = bilby.core.prior.Uniform(name='cos_iota',minimum=-1,maximum=1)
PriorsX['psi']      = bilby.core.prior.Uniform(name='psi',minimum=0,maximum=np.pi)
PriorsX['t_coal']   = bilby.core.prior.Uniform(name='t_coal',minimum=-1.,maximum=1.)
PriorsX['phi_coal'] = bilby.core.prior.Uniform(name='phi_coal',minimum=-np.pi,maximum=np.pi)

PriorsX['sx1']      = bilby.core.prior.Uniform(name='sx1',minimum=-0.98,maximum=0.98)
PriorsX['sy1']      = bilby.core.prior.Uniform(name='sy1',minimum=-0.98,maximum=0.98)
PriorsX['sz1']      = bilby.core.prior.Uniform(name='sz1',minimum=-0.98,maximum=0.98)

PriorsX['sx2']      = bilby.core.prior.Uniform(name='sx2',minimum=-0.98,maximum=0.98)
PriorsX['sy2']      = bilby.core.prior.Uniform(name='sy2',minimum=-0.98,maximum=0.98)
PriorsX['sz2']      = bilby.core.prior.Uniform(name='sz2',minimum=-0.98,maximum=0.98)

PriorsX["S1"]		= bilby.core.prior.Uniform(name='S1',minimum=0.,maximum=0.98)
PriorsX["phi1"]		= bilby.core.prior.Uniform(name='phi1',minimum=0.,maximum=2*np.pi)
PriorsX['theta1']   = bilby.core.prior.Sine(name='theta1',minimum=0,maximum=np.pi)

PriorsX["S2"]		= bilby.core.prior.Uniform(name='S2',minimum=0.,maximum=0.98)
PriorsX["phi2"]		= bilby.core.prior.Uniform(name='phi2',minimum=0.,maximum=2*np.pi)
PriorsX['theta2']   = bilby.core.prior.Sine(name='theta2',minimum=0,maximum=np.pi)

PriorsX['chi_s']  = bilby.core.prior.Uniform(name='chi_s',minimum=-0.98,maximum=0.98)
PriorsX['chi_a']  = bilby.core.prior.Uniform(name='chi_a',minimum=-0.98,maximum=0.98)

def Check_Priors(FreeParams,plot_name=None,new_priors=None,plot=False):
	
	priors = bilby.core.prior.PriorDict()#conversion_function=MassConstraint)
	for fp in FreeParams:
		priors[fp] = PriorsX[fp] 
	if(new_priors!=None):
		for key in new_priors.keys():
			x_vec, y_vec = new_priors[key]
			priors[key] = bilby.core.prior.Interped(name=key,xx=x_vec,yy=y_vec,minimum=min(x_vec),maximum=max(x_vec))

	if plot: 
		print("Plotting Priors ...")

		fig = plt.figure(figsize=(10,8))
		ndim = len(FreeParams)
		N = np.sqrt(ndim)
		if N%int(N)>0: N = int(N)+1
		else: N = int(N)
		for i in range(ndim):
			plt.subplot(N,N,i+1)
			key = FreeParams[i]
			x = np.linspace(priors[key].minimum , priors[key].maximum, 1000 )
			y = priors[key].prob(x)
			y = y/max(y)

			idxs = np.where(y>0)[0]
			xx = x[idxs]

			zeros = np.ones(len(x))*min(y)

			plt.plot(x,y,'k-',lw=2)
			plt.fill_between(x,zeros,y,color='grey',alpha=0.5)
			plt.grid(alpha=0.3)
			plt.xlabel(key)
		plt.tight_layout()
		if plot_name==None: plot_name = "priors.png"
		fig.savefig(plot_name)
		print("plot of priors OK!")

	return priors

def get_samples(Data,Tensors,detectors,gwprms,FreeParams,
				approx,Method,sampler,npoints,nwalkers,ntemps,nburn,
				verbose,new_priors,limits=None,pos0=None,npool=None,
				remove_out=True,output_name=None,save_bilby_path=True,
				enable_jax_waveforms=True,bilby_path="outputs_bilby/",**kwargs):
	priors = {}
	for fp in FreeParams:
		priors[fp] = PriorsX[fp]
	if(new_priors!=None):
		for key in new_priors.keys():
			x_vec, y_vec = new_priors[key]
			priors[key] = bilby.core.prior.Interped(name=key,xx=x_vec,yy=y_vec,minimum=min(x_vec),maximum=max(x_vec))

	theta_keys, gw_values = zip(*gwprms.items())
	if enable_jax_waveforms:
		GwSignal = wf.load_waveforms(theta_keys,approx,**kwargs)[0]
		Gw_Signal = lambda *args: GwSignal[approx](*args)
	else:
		h_func = wf.build_waveform_strain_lal(theta_keys,approx,**kwargs)
		Gw_Signal = lambda *args: jnp.asarray(h_func(*args, approx,**kwargs))

	if(Method=="Exact"): likelihood = Exact_likelihood([Data,detectors,gwprms,Gw_Signal,FreeParams])
	else: likelihood = DALI_likelihood([gwprms,Tensors,FreeParams,Method])

	if(sampler=="grid"):
		return grd.Call_Grid(likelihood,priors,npoints,limits)
	else:
		idx = np.random.uniform(0,int(1.e5))
		os.makedirs(bilby_path,exist_ok=True)
		if output_name == None:
			output_name = f"{len(FreeParams)}_out_{idx}"
		outdir = f"{bilby_path}{output_name}/"
		print("\n\n Running Sampler (%s)...\n" % sampler)

		if pos0 is not None:
			print(f"\nsampler: {sampler}")
			print("FreeParams:",FreeParams)
			print(f"ndim: {len(FreeParams)}")
			print(f"nsteps: {npoints}")
			print(f"nwalkers: {nwalkers}")
			print(f"ntemps: {ntemps}")
			print(f"npool: {npool}")
			print("shape(pos0):",np.shape(pos0),'\n')

			res = bilby.run_sampler(likelihood,priors,sampler=sampler,npoints=npoints,nwalkers=nwalkers,\
									nsteps=npoints,nburn=nburn,outdir=outdir,npool=npool, verbose=verbose,
									store_walkers=False,resume=True,save=save_bilby_path,check_point_plot=False,ntemps=ntemps,pos0=pos0,**kwargs)
		else:
			res = bilby.run_sampler(likelihood,priors,sampler=sampler,npoints=npoints,nwalkers=nwalkers,\
									nsteps=npoints,nburn=nburn,outdir=outdir,npool=npool, verbose=verbose,
									store_walkers=False,resume=True,save=save_bilby_path,check_point_plot=False,ntemps=ntemps,**kwargs)

		if(os.path.isdir(outdir) and remove_out): os.system('rm -R %s' % outdir)

		Samples = []
		for fp in FreeParams:
			Samples.append( res.posterior[fp] )
		try:
			Evidence = [res.log_evidence, res.log_evidence_err]
		except:
			Evidence = [None,None]

		likelihood_times = np.array( likelihood.get_times() )
		return [np.transpose(Samples),priors,likelihood_times,Evidence]