from jax import jit, lax, config, debug
import jax
import jax.numpy as jnp
import numpy as np
import GWDALI.lib.AngTransf as geo
import GWDALI.lib.TaylorF2 as Tf2
import GWDALI.lib.IMRPhenomA as phnA
import GWDALI.lib.IMRPhenomB as phnB
import GWDALI.lib.IMRPhenomC as phnC
import GWDALI.lib.IMRPhenomD as phnD
import GWDALI.lib.IMRPhenomHM as phnHM

import warnings
warnings.filterwarnings("ignore")

'''
# Approximants:
	1. TaylorF2_Spinless (3.5 PN  - No Spins)
	2. TaylorF2_Spinless_0PN (0 PN - No Spins)
	3. TaylorF2 (3.5 PN - Spins)
	4. TaylorF2_ISCO (3.5 PN - Merger Cutoff)
	5. IMRPhenomA
	6. IMRPhenomB
	7. IMRPhenomC
	8. IMRPhenomD
	9. IMRPhenomHM
'''

PI = 3.141592653589793238462643383279502884
gE = 0.577215664901532860606512090082402431
Msun = 1.988409870698050731911960804878414216e30
pc = 3.085677581491367278913937957796471611e16
G = 6.67430e-11
c = 299792458.

pi2 = PI*PI

rad = PI/180
w_earth = 2*PI/(24.*3600)

R_earth_equatiorial = 6.3781e6
R_earth_polar = 6.3568e6
R_earth = 6.378e6 # meters

A0 = 7.806521525937888e-23
# A0 := sqrt(5/24)*(G*Msun/c^3)^(5./6) / (pi^(2/3) * Gpc / c)

GMc2 = 1476.6250380501247 # GM/c^2
GMc3 = 4.925490947641267e-06 # GM/c^3

#===============================================================

@jit
def PatternFunc(alpha,beta,psi,shape): # L-shape Detector
	u = jnp.cos(beta)
	fp = 0.5*(1+u*u)*jnp.cos(2*alpha)
	fx = -u*jnp.sin(2*alpha)
	Fp = fp*jnp.cos(2*psi) + fx*jnp.sin(2*psi)
	Fx =-fp*jnp.sin(2*psi) + fx*jnp.cos(2*psi)
	Omega = shape*PI/180
	return Fp*jnp.sin(Omega), Fx*jnp.sin(Omega)

@jit
def PN_time_corrections(m1,m2,freq):
	# Equation (3.8b) of arXiv:0907.0700 
	# t(f) up to 3.5 PN 
	# spinless!

	M = m1+m2
	eta = m1*m2/M**2

	v  = 2*PI*M*GMc3*freq
	v2 = v*v
	v3 = v2*v
	v4 = v2*v2
	v5 = v3*v2
	v6 = v3*v3
	v7 = v4*v3
	v8 = v4*v4

	eta2 = eta*eta
	eta3 = eta2*eta

	a0 = 1.
	a1 = 0.
	a2 = 743./252 + 11*eta/3
	a3 = -32.*PI/5
	a4 = 3058673./508032 + 5429.*eta/504 + 617./72
	a5 = -(7728./252 - 13.*eta/3)*PI
	a6 = -10052469856691./23471078400 + 128./3*PI**2 + 6848.*jnp.euler_gamma/105 + \
	(3147553127./3048192 - 451.*PI**2/12)*eta - 15211./1728*eta2 + 25565.1296*eta3 + 3424.*jnp.log(16*v2)/105
	a7 = (-15419335./127008 - 75703.*eta/756 + 14809*eta2/378)*PI

	Corrections = a0 + a1*v + a2*v2 + a3*v3 + a4*v4 + a5*v5 + a6*v6 + a7*v7 # 3.5 PN
	Coeff =  - 5*M*GMc3 / (256*eta*v8) # Leading Order (0PN)

	return Coeff*Corrections

@jit
def get_FpFx(iota,psi,RA,Dec,det_conf,Mc,eta,freq):
	alpha0 = RA*rad
	beta0 = Dec*rad
	lon, lat, rot, shape = det_conf

	alpha_det, beta_det, psi_det = geo.AngTransf(iota,psi,RA,Dec,lon,lat,rot)
	Fp, Fx = PatternFunc(alpha_det, beta_det, psi_det, shape)
	return Fp, Fx

@jit
def get_FpFx_ER(iota,psi,RA,Dec,det_conf,Mc,eta,freq):
	alpha0 = RA*rad
	beta0 = Dec*rad
	lon, lat, rot, shape = det_conf

	M = Mc/eta**(3./5)
	dM = jnp.sqrt(1-4*eta)
	m1 = .5*M*(1.+dM)
	m2 = .5*M*(1.-dM)

	t_shift = PN_time_corrections(m1,m2,freq) # (t_shift<=0) 3.5PN time(f) [arXiv:0907.0700 ]
	RA_mod = RA - w_earth*t_shift/rad # Effect of Earth Rotation
	
	alpha_det, beta_det, psi_det = geo.AngTransf(iota,psi,RA_mod,Dec,lon,lat,rot)
	Fp, Fx = PatternFunc(alpha_det, beta_det, psi_det, shape)
	return Fp, Fx

hphx_TaylorF2              = Tf2.hphx_TaylorF2
hphx_TaylorF2_ISCO         = Tf2.hphx_TaylorF2_ISCO
hphx_TaylorF2_Spinless     = Tf2.hphx_TaylorF2_Spinless
hphx_TaylorF2_Spinless_0PN = Tf2.hphx_TaylorF2_Spinless_0PN

hphx_IMRPhenomA  = phnA.hphx_IMRPhenomA
hphx_IMRPhenomB  = phnB.hphx_IMRPhenomB
hphx_IMRPhenomC  = phnC.hphx_IMRPhenomC
hphx_IMRPhenomD  = phnD.hphx_IMRPhenomD
hphx_IMRPhenomHM = phnHM.hphx_IMRPhenomHM

APPROXIMANTS = [
	'TaylorF2_Spinless_0PN', # 0 PN  + cutoff in f=f_ISCO
	'TaylorF2_Spinless', # 3.5 PN  + cutoff in f=f_ISCO
	'TaylorF2', # 3.5 PN (no frequency cutoff)
	'TaylorF2_ISCO', # 3.5 PN + cutoff in f=f_ISCO
	'IMRPhenomA',
	'IMRPhenomB',
	'IMRPhenomC',
	'IMRPhenomD',
	'IMRPhenomHM',
]

@jit
def get_time_delay(det_conf_a,det_conf_b,iota,psi,RA,Dec):
	beta_a = geo.AngTransf(iota,psi,RA,Dec,*det_conf_a[:-1])[1]
	beta_b = geo.AngTransf(iota,psi,RA,Dec,*det_conf_b[:-1])[1]
	
	return R_earth/c * (jnp.cos(beta_b) - jnp.cos(beta_a))

#==============================================================================================
# LAL WAVEFORMS
#==============================================================================================
Gpc = pc*1.e9
wrn =\
'''
=======================================================
    WARNING: lalsuite or lalsimulation not instaled!
    Try to install them with: 
    \033[1m conda install -c conda-forge lalsuite \033[0m
    \033[1m conda install -c conda-forge lalsimulation \033[0m
======================================================='''
try: 
    import lal
    import lalsimulation as lalsim
except:
    print(wrn)
    
from scipy.interpolate import interp1d

def interp_complex(freq0, h, freq):
	cond = True
	if cond:
		amp = np.abs(h)
		#amp  = np.sqrt(np.real(h*np.conj(h)))
		phase = np.unwrap(np.angle(h))
		amp_i = interp1d(freq0, amp, bounds_error=False, fill_value=0.0)
		phase_i = interp1d(freq0, phase, bounds_error=False, fill_value=0.0)
		h = amp_i(freq) * np.exp(1j * phase_i(freq))
	else:
		hr = np.real(h); funcR = interp1d(freq0,hr,bounds_error=False,fill_value=0.) ; hr = funcR(freq)
		hi = np.imag(h); funcI = interp1d(freq0,hi,bounds_error=False,fill_value=0.) ; hi = funcI(freq)
		h = hr + 1.j*hi
	return h

def hphx_lal(dL, iota, phi0, t0, Mc, eta, sx1, sy1, sz1, sx2, sy2, sz2, freq, approx, **kwargs):
	approx_lal = lalsim.GetApproximantFromString(approx)
	
	deltaM = np.sqrt(1.-4*eta)
	M = Mc/eta**(3./5)
	m1 = 0.5*M * (1. + deltaM) * Msun
	m2 = 0.5*M * (1. - deltaM) * Msun

	f = np.array(freq)
	f_low, f_max = np.min(f), np.max(f)
	
	dF = 0.125
	if "dF" in kwargs.keys():
		if kwargs["dF"] != None:
			dF = kwargs["dF"]

	dL_SI = dL * Gpc # Gpc to meters

	phi_ref = 0.
	fRef = 1.0 if approx in ["IMRPhenomD","IMRPhenomHM"] else 0.0

	hp,hx = lalsim.SimInspiralChooseFDWaveform(
		    m1, m2,
		    sx1, sy1, sz1, 
		    sx2, sy2, sz2,
		    dL_SI,
		    iota, phi_ref, 0., 0., 0.,
		    dF, f_low, f_max, fRef,
		    None,
		    approx_lal
		)
	hp = hp.data.data
	hx = hx.data.data

	N = len(hp)
	freq0 = np.linspace(0,dF*(N-1),N)

	Exp0 = np.exp( -1.j * (2*PI*freq0*t0 - phi0) )

	hp_new = interp_complex(freq0, hp*Exp0, f) #* Exp0
	hx_new = interp_complex(freq0, hx*Exp0, f) #* Exp0
	return hp_new, hx_new

# ===============================================================
# ======================= BUILD WAVEFORM =======================
# ===============================================================
from GWDALI.lib.ParamsTransform import build_transform

def build_waveform_strain(theta_keys,approx,**kwargs):
	transform = build_transform(theta_keys)

	earth_rotation = kwargs.get("ER", False)
	FpFx_func = ( get_FpFx_ER if earth_rotation else get_FpFx )

	@jit
	def Gw_Signal(theta, det_conf, freq):
		prms = transform(theta)
		dL,iota,psi,phi_coal,RA,Dec,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2 = prms
		args_hphx = [dL,iota,phi_coal,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq]
		args_FpFx = [RA,Dec,psi]

		hp, hx = globals()[f"hphx_{approx}"](*args_hphx)

		det_conf_a, det_conf_ref = det_conf
		tau_ab = get_time_delay(det_conf_a,det_conf_ref,iota,psi,RA,Dec)

		Fp, Fx = FpFx_func(iota,psi,RA,Dec,det_conf_a,Mc,eta,freq)
		return (Fp*hp + Fx*hx) * jnp.exp(2.j * jnp.pi * freq * tau_ab)
	return Gw_Signal
# ===============================================================

def build_waveform_hphx(theta_keys,approx):
    transform = build_transform(theta_keys)
    @jit
    def Gw_hphx(theta, freq):
        prms = transform(theta)
        dL,iota,psi,phi_coal,RA,Dec,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2 = prms
        args_hphx = [dL,iota,phi_coal,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq]
        return globals()[f"hphx_{approx}"](*args_hphx)
    return Gw_hphx
# ===============================================================

def load_waveforms(theta_keys,approx,**kwargs):
	Gw_hphx = {}
	Gw_Signal = {}
	Gw_Signal_Real = {}
	Gw_Signal_Imag = {}

	gw_strain = build_waveform_strain(theta_keys,approx,**kwargs)
	Gw_hphx[approx]        = build_waveform_hphx(theta_keys,approx)   
	Gw_Signal[approx]      = gw_strain 
	Gw_Signal_Real[approx] = lambda *args: jnp.real(gw_strain(args,args[-2],args[-1]))
	Gw_Signal_Imag[approx] = lambda *args: jnp.imag(gw_strain(args,args[-2],args[-1]))

	return Gw_Signal, Gw_hphx, Gw_Signal_Real, Gw_Signal_Imag

#========================================================

def build_waveform_hphx_lal(theta_keys,approx,**kwargs):
	transform = build_transform(theta_keys)
	def Gw_hphx(theta, freq, approx, **kwargs):
		prms = np.array( transform(theta) ,dtype=float )
		dL,iota,psi,phi_coal,RA,Dec,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2 = prms
		args_hphx = [dL,iota,phi_coal,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq]
		return hphx_lal(*args_hphx,approx,**kwargs)
	return Gw_hphx

def build_waveform_strain_lal(theta_keys,approx,**kwargs):
	transform = build_transform(theta_keys)
	
	earth_rotation = kwargs.get("ER", False)
	FpFx_func = ( get_FpFx_ER if earth_rotation else get_FpFx )

	def Gw_Signal(theta, det_conf, freq, approx, **kwargs):
		prms = np.array( transform(theta) ,dtype=float )
		dL,iota,psi,phi_coal,RA,Dec,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2 = prms
		hphx_prms = [dL, iota, phi_coal, t_coal, Mc, eta, sx1, sy1, sz1, sx2, sy2, sz2]
		hp, hx = hphx_lal(*hphx_prms, freq, approx, **kwargs)

		det_conf_a, det_conf_ref = det_conf
		tau_ab = get_time_delay(det_conf_a,det_conf_ref,iota,psi,RA,Dec)

		Fp, Fx = FpFx_func(iota,psi,RA,Dec,det_conf_a,Mc,eta,freq)
		return (Fp*hp + Fx*hx) * jnp.exp(2.j * jnp.pi * freq * tau_ab)
	return Gw_Signal