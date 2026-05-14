from jax import jit
import jax.numpy as jnp
import numpy as np

PI = 3.141592653589793238462643383279502884
gE = 0.577215664901532860606512090082402431
Msun = 1.988409870698050731911960804878414216e30
pc = 3.085677581491367278913937957796471611e16
G = 6.67430e-11
c = 299792458.
Mpc = pc*1.e6

GMc3 = G*Msun/c**3 #GMc3 = G*Msun/c**3

A0 = 7.806521525937888e-23

@jit # from 1PN to 3.5 PN
def PhaseTf2(Mc,eta,freq):
	M = Mc/eta**(3./5)
	eta2 = eta*eta
	eta3 = eta2*eta
	pi2 = PI*PI
	Mf = GMc3*M*freq
	piMf = PI*Mf

	#  .5 PN (a1) ; 1 PN (a2) ; 1.5 PN (a3) ; 2 PN (a4) ; 2.5 PN (a5) ; 3 PN (a6) ; 3.5 PN (a7) 
	a1 = 0. 
	a2 = 3715./756 + 55./9*eta 
	a3 = -16*PI 
	a4 = 15293365./508032 + 27145./504*eta + 3085./72*eta2 
	a5 = PI * ( 38645./756 - 65/9*eta ) * ( 1. + jnp.log( 6**(3./2)*PI*Mf) ) 
	a6 = 11583231236531./4694215680 - 640./3*pi2 - 6848./21*gE\
		 + ( -15737765635./3048192 + 2255./12*pi2 )*eta \
		 + 76055./1728*eta2 - 127825./1296*eta3\
		 - 6848./63*jnp.log( 64*PI*Mf ) 
	a7 = PI*( 77096675./254016 + 378515./1512*eta - 74045./756*eta2 )
	
	As = [0,a1,a2,a3,a4,a5,a6,a7]
	v = piMf**(1./3)

	Phi = 0.0 ; pfN = 3./(128*eta*v**5)
	for k in range(8): Phi += As[k]*v**k
	
	return pfN*Phi

@jit # from 1PN to 3.5 PN + Spins
def PhaseTf2_Spins(Mc,eta,chi1,chi2,freq):
	M = Mc/eta**(3./5)
	eta2 = eta*eta
	eta3 = eta2*eta
	pi2 = PI*PI

	m1 = 0.5*(Mc/eta**(3./5)) * (1 + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1 - jnp.sqrt(1-4*eta))
	m1M = m1/M ; m2M = m2/M
	deltaM = m1M - m2M
	chi_a = (chi1-chi2)/2
	chi_s = (chi1+chi2)/2

	Mf = GMc3*M*freq
	piMf = PI*Mf
	v = piMf**(1./3)

	varphi = [0]*8

	#varphi[1] += 1;
	varphi[2] += 5.*(74.3/8.4 + 11.*eta)/9.;
	varphi[3] += -16.*PI;
	varphi[4] += 5.*(3058.673/7.056 + 5429./7.*eta+617.*eta*eta)/72.;
	varphi[5] += 5./9.*(772.9/8.4-13.*eta)*PI;
	#vlogv[5] += 5./3.*(772.9/8.4-13.*eta)*PI;
	varphi[6]     =  11583.231236531/4.694215680 - 640./3.*pi2 - 684.8/2.1*gE + eta*(-15737.765635/3.048192 + 225.5/1.2*pi2) + eta2*76.055/1.728 - eta3*127.825/1.296;
	#vlogvarphi[6] += -684.8/2.1;
	varphi[7] += PI*(770.96675/2.54016 + 378.515/1.512*eta - 740.45/7.56*eta*eta);
	#======================================================================
	# Spinning Corrections to varphi (from lalsimulation):
	#======================================================================
	varphi[7] += ( m1M*(-17097.8035/4.8384+eta*28764.25/6.72+eta2*47.35/1.44 + m1M*(-7189.233785/1.524096+eta*458.555/3.024-eta2*534.5/7.2)) )*chi1 \
	+ ( m2M*(-17097.8035/4.8384+eta*28764.25/6.72+eta2*47.35/1.44 + m2M*(-7189.233785/1.524096+eta*458.555/3.024-eta2*534.5/7.2)) )*chi2

	varphi[6] += ( PI*m1M*(1490./3. + m1M*260.) )*chi1 \
	+ ( PI*m2M*(1490./3. + m2M*260.) )*chi2 \
	+ ( (326.75/1.12 + 557.5/1.8*eta)*eta )*chi1*chi2 \
	+ ( ( (4703.5/8.4+2935./6.*m1M-120.*m1M**2)*m1M**2 ) + ( (-4108.25/6.72-108.5/1.2*m1M+125.5/3.6*m1M**2)*m1M**2 ) )*chi1*chi1 \
	+ ( ( (4703.5/8.4+2935./6.*m2M-120.*m2M**2)*m2M**2 ) + ( (-4108.25/6.72-108.5/1.2*m2M+125.5/3.6*m2M**2)*m2M**2 ))*chi2*chi2

	varphi[5] += ( -m1M*(1391.5/8.4-m1M*(1.-m1M)*10./3.+ m1M*(1276./8.1+m1M*(1.-m1M)*170./9.)) )*chi1 \
	+  ( -m2M*(1391.5/8.4-m2M*(1.-m2M)*10./3.+ m2M*(1276./8.1+m2M*(1.-m2M)*170./9.)) )*chi2

	varphi[4] += ( 247./4.8*eta )*chi1*chi2 \
	+ ( -721./4.8*eta )*chi1*chi2 \
	+ ( ( -720./9.6*m1M**2 ) + ( 1./9.6*m1M**2 ) )*chi1*chi1 \
	+ ( ( -720./9.6*m2M**2 ) + ( 1./9.6*m2M**2 ) )*chi2*chi2 \
	+ ( ( 240./9.6*m1M**2 ) + ( -7./9.6*m1M**2 ) )*chi1*chi1 \
	+ ( ( 240./9.6*m2M**2 ) + ( -7./9.6*m2M**2 ) )*chi2*chi2

	varphi[3] += ( m1M*(25.+38./3.*m1M) )*chi1 +  ( m2M*(25.+38./3.*m2M) )*chi2
	#======================================================================

	VarPhi = [x for x in varphi]
	VarPhi[5] *= (1.+jnp.log(piMf))
	#VarPhi[5] *= (1.+jnp.log(piMf * 6**1.5)) # Why do we have to include 6**1.5
	VarPhi[6] -= 6848*jnp.log(64*piMf)/63

	pfN = 3./(128*eta*v**5)
	Phase = v*0
	for k in range(8): Phase += VarPhi[k]*v**k
	Phase *= pfN
	return Phase

@jit # Phase TaylorF2 (0PN)
def PhaseTf2_0PN(phi0,t0,Mc,eta,freq):
	M = Mc/eta**(3./5)

	f_isco = 1./(6*jnp.sqrt(6)*M) / GMc3 / jnp.pi
	cutoff = freq<f_isco
	Phase = 2*PI*freq*t0 - phi0 - PI/4 + (3./128)*(PI*GMc3*Mc*freq)**(-5./3)
	return Phase

@jit # Amplitude TaylorF2
def Amp_Tf2(dL,Mc,eta,freq):

	M = Mc/eta**(3./5)

	h0 = A0 * Mc**(5./6) / dL
	#f_isco = 1./(6*jnp.sqrt(6)*M) / GMc3 / jnp.pi
	#cutoff = freq<f_isco
	return h0*freq**(-7./6)

@jit # TaylorF2 (0PN) Waveform
def hphx_TaylorF2_Spinless_0PN(dL,iota,phi0,t0,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq):

	M = Mc/eta**(3./5)
	f_isco = 1./(6*jnp.sqrt(6)*M) / GMc3 / jnp.pi
	cutoff = freq<f_isco

	Phase = PhaseTf2_0PN(phi0,t0,Mc,eta,freq)
	Amp = Amp_Tf2(dL,Mc,eta,freq)*cutoff

	h0 = Amp*jnp.exp(1.j*Phase)

	u = jnp.cos(iota)
	gp = 0.5*(1+u*u)
	gx = -1.j*u

	hp = gp*h0
	hx = gx*h0

	return hp, hx

@jit # TaylorF2 (3.5 PN) Waveform
def hphx_TaylorF2_Spinless(dL,iota,phi0,t0,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq):
	
	Phase_corr = PhaseTf2(Mc,eta,freq)
	Phase = PhaseTf2_0PN(phi0,t0,Mc,eta,freq) + Phase_corr
	
	M = Mc/eta**(3./5)
	f_isco = 1./(6*jnp.sqrt(6)*M) / GMc3 / jnp.pi
	cutoff = freq<f_isco

	Amp = Amp_Tf2(dL,Mc,eta,freq)*cutoff
	h0 = Amp*jnp.exp(1.j*Phase)

	u = jnp.cos(iota)
	gp = 0.5*(1+u*u)
	gx = -1.j*u

	hp = gp*h0
	hx = gx*h0

	return hp, hx

@jit # TaylorF2 (3.5 PN) Waveform + Spins
def hphx_TaylorF2(dL,iota,phi0,t0,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq):
	
	Phase_spins = PhaseTf2_Spins(Mc,eta,sz1,sz2,freq)
	Phase_0PN = 2*PI*freq*t0 - phi0 - PI/4 + (3./128)*(PI*GMc3*Mc*freq)**(-5./3)

	Phase = -( Phase_0PN + Phase_spins ) + PI
	
	M = Mc/eta**(3./5)
	f_isco = 1./(6*jnp.sqrt(6)*M) / GMc3 / jnp.pi
	cutoff = freq<f_isco

	Amp = Amp_Tf2(dL,Mc,eta,freq)#*cutoff
	h0 = Amp*jnp.exp(1.j*Phase)

	u = jnp.cos(iota)
	gp = 0.5*(1+u*u)
	gx = -1.j*u

	hp = gp*h0
	hx = gx*h0

	return hp, hx

@jit # TaylorF2 (3.5 PN) Waveform + Spins + Cutoff
def hphx_TaylorF2_ISCO(dL,iota,phi0,t0,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq):
	
	M = Mc/eta**(3./5)
	f_isco = 1./(6*jnp.sqrt(6)*M) / GMc3 / jnp.pi
	cutoff = freq<f_isco

	hp, hx = hphx_TaylorF2(dL,iota,phi0,t0,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq)

	return hp*cutoff, hx*cutoff
