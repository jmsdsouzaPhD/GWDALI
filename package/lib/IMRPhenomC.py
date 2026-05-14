#=========================================================#
# IMRPhenomC [arXiv:1005.3306v3]
#=========================================================#
# Check LALSimIMRPhenomC.c (line 303)
# https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_i_m_r_phenom_c_8c_source.html#l00303
# Check LALSimIMRPhenomC_internals.c (line 373)
# https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_i_m_r_phenom_c__internals_8c_source.html#l00373

from jax import jit
import jax.numpy as jnp
import numpy as np
#'from jax.scipy.interpolate import interp1d

PI = 3.141592653589793238462643383279502884
gE = 0.577215664901532860606512090082402431
Msun = 1.988409870698050731911960804878414216e30
pc = 3.085677581491367278913937957796471611e16
G = 6.67430e-11
c = 299792458.
Mpc = pc*1.e6

GMc3 = G*Msun/c**3
A0 = 7.806521525937888e-23

pi2 = PI*PI

# Constants to estimate final spin
s4 = -0.129
s5 = -0.384
tt0 = -2.686
t2 = -3.454
t3 = 2.353

#------------------------------------------#------------------------------------------
# From Table II of PhysRevD.82.064016:
#------------------------------------------#------------------------------------------
# Coeffs for alpha_k

Zeta_a1 = [-2.417e-3, -1.093e-3, -1.917e-2, 7.267e-2, -2.504e-1]
Zeta_a2 = [5.962e-1, -5.600e-2, 1.520e-1, -2.970e0, 1.312e1]
Zeta_a3 = [-3.283e1, 8.859e0, 2.931e1, 7.954e1, -4.349e2]
Zeta_a4 = [1.619e2, -4.702e1, -1.751e2, -3.225e2, 1.587e3]
Zeta_a5 = [-6.320e2, 2.463e2, 1.048e3, 3.355e2, -5.115e3]
Zeta_a6 = [-4.809e1, -3.643e2, -5.215e2, 1.870e3, 7.354e2]
Zeta_g1 = [4.149e0, -4.070e0, -8.752e1, -4.897e1, 6.665e2]
Zeta_d1 = [-5.472e-2, 2.094e-2, 3.554e-1, 1.151e-1, 9.640e-1]
Zeta_d2 = [-1.235e0, 3.423e-1, 6.062e0, 5.949e0, -1.069e1]

#------------------------------------------#------------------------------------------

@jit
def Lorentzian(f,f0,sig):
	sig2 = sig*sig
	df = (f-f0)
	df2 = df*df
	return sig2/( df2 + sig2/4 )

@jit
def func_Q(a):
	q1, q2, q3 = 0.7, 1.4187, -0.4990
	return q1 + q2*(1-a)**q3

@jit
def func_fRD(a): # Dimensionless Ringdown frequency
	k1, k2, k3 = 1.5251, -1.1568, 0.1292
	return (k1 + k2*(1-a)**k3)/(2*PI)

#================================================================
# Transition functions:
#================================================================
@jit
def wMinus(Mf,d,f0):
	df = (Mf-f0)
	w = 0.5*( 1. - jnp.tanh(4*df/d) )
	return w

@jit
def wPlus(Mf,d,f0):
	df = (Mf-f0)
	w = 0.5*( 1. + jnp.tanh(4*df/d) )
	return w

@jit
def diff_wMinus(Mf,d,f0):
	df = (Mf-f0)
	u = jnp.tanh(4*df/d)
	w = -(4/d)*( 1. - u*u )
	return w

@jit
def diff_wPlus(Mf,d,f0):
	df = (Mf-f0)
	u = jnp.tanh(4*df/d)
	w = (4/d)*( 1. - u*u )
	return w

#================================================================
# Phase Functions: Inspiral-Merger-RingDown
#================================================================
@jit
def PhaseC_SPA(eta,Alpha,v): # Stationary Phase Apprximation (SPA)
	Summation = v*0 ; v5 = v**5
	piMf = v3 = v*v*v
	Alpha[5] *= (1.+jnp.log(piMf))
	Alpha[6] -= 6848*jnp.log(64*piMf)/63
	for k in range(8):
		Summation += Alpha[k]*v**k
	N = 3./( 128.*eta*v**5 )
	return N*Summation - PI/4

@jit
def PhaseC_PM(eta,Alpha,f): # Pre-Merger Phase
	al1, al2, al3, al4, al5, al6 = Alpha
	return ( al1*f**(-5./3) + al2/f + al3*f**(-1./3) + al4 + al5*f**(2./3) + al6*f )/eta

@jit
def PhaseC_RD(f,Beta):
	beta1, beta2 = Beta
	return beta1 + beta2*f

#================================================================
# Derivatives of Phase Functions: Inspiral-Merger-RingDown
#================================================================
@jit
def diff_PhaseC_SPA(eta,Alpha,v): # Stationary Phase Apprximation (SPA)
	Summation = v*0 ; v5 = v**5
	piMf = v3 = v*v*v
	diff_Alpha = [x*0 for x in Alpha]
	diff_Alpha[5] += PI/v3
	diff_Alpha[6] -= 6848*64*PI/(63*v3)

	Alpha[5] *= (1.+jnp.log(piMf))
	Alpha[6] -= 6848*jnp.log(64*piMf)/63
	for k in range(8):
		Summation += ((k-5)*Alpha[k]*v**(k-1) + diff_Alpha[k]*v**k)
	N = 3./( 128.*eta*v5 )
	return N*Summation

@jit
def diff_PhaseC_PM(eta,Alpha,f): # Derivative of Pre-Merger Phase
	al1, al2, al3, al4, al5, al6 = Alpha
	return ( (-5./3)*al1*f**(-5./3) - al2/f + (-1./3)*al3*f**(-1./3) + (2./3)*al5*f**(2./3) + al6*f )/eta/f

@jit
def diff_PhaseC_RD(f,Beta):
	return Beta[1]

#================================================================

@jit
def diff_PhaseC_PM(eta,Alpha,f): # Derivative of Pre-Merger Phase
	al1, al2, al3, al4, al5, al6 = Alpha
	return ( (-5./3)*al1*f**(-5./3) - al2/f + (-1./3)*al3*f**(-1./3) + (2./3)*al5*f**(2./3) + al6*f )/eta/f

@jit
def Phase_IMRPhenomC(t0,phi0,Mc,eta,chi1,chi2,freq):

	m1 = 0.5*(Mc/eta**(3./5)) * (1 + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1 - jnp.sqrt(1-4*eta))

	chi2 = chi1
	M    = m1+m2

	eta2 = eta*eta
	eta3 = eta2*eta

	GM   = M*GMc3
	f    = GM*freq
	piMf = PI*f
	
	xi = (m1*chi1 + m2*chi2)/M
	chi_s = chi1 + chi2
	chi_p = chi1*chi2
	xi2 = xi*xi
	xi3 = xi2*xi

	#=======================================
	# Stationary Phase Apprximation (SPA)
	#=======================================
	aIns = [ freq*0 for i in range(8)]

	aIns[0] += 1.
	aIns[1] += 0.
	aIns[2] += 3715./756. + 55.*eta/9
	aIns[3] += -16*PI + 113*xi/3 - (38*eta/3)*chi_s
	aIns[4] += 15293365./508032 - 50*xi2 + eta*(27145/504 + 5*chi_p/4) + 3085*eta2/72
	aIns[5] += ( PI*(38645/756 - 65*eta/9) - xi*(735505./2268 + 130*eta/9) \
				 + chi_s*(12850*eta/81 + 170*eta2/9) - 10*xi3/3 + 10*eta*xi*chi_p)
	aIns[6] += 11583231236531./4694215680 - 640*pi2/3 - 6848*gE/21 \
				 + eta*(2255*pi2/12 - 15737765635./3048192) + 76055*eta2/1728 \
				 - 127825*eta3/1296 + 2920*PI*xi/3 - (175.-1490*eta)*xi2/3 \
				 - (1120*PI/3 - 1085*xi/3)*eta*chi_s + (26945*eta/336 - 2365*eta2/6)*chi_p
	aIns[7] += PI*(77096675./254016. + 378515.*eta/1512. - 74045.*eta2/756.) \
				 - xi*(20373952415./3048192. + 150935.*eta/224. - 578695.*eta2/432)\
				 + chi_s*( 4862041225.*eta/1524096. + 1189775.*eta2/1008. -71705.*eta3/216\
				 - 830.*eta*xi2/3 + 35*eta2*chi_p/3) - 560*PI*xi2\
				 + 20*PI*eta*chi_p + xi3*(94555./168 - 85*eta) + xi*chi_p*( 39665.*eta/168 + 255*eta2 )
	
	Summation = freq*0 ; v = piMf**(1./3)
	for k in range(8):
		Summation += aIns[k]*v**k

	N = 3./( 128.*eta*v**5 )
	Phase_SPA = N*Summation - PI/4 
	
	#=======================================
	# Pre-Merger (PM)
	#=======================================

	eta_xi = [xi,xi2,eta*xi,eta,eta2]
	aPM = [0]*6
	aPM[0] = jnp.sum( jnp.array( [ Zeta_a1[k]*eta_xi[k] for k in range(5) ] ) )
	aPM[1] = jnp.sum( jnp.array( [ Zeta_a2[k]*eta_xi[k] for k in range(5) ] ) )
	aPM[2] = jnp.sum( jnp.array( [ Zeta_a3[k]*eta_xi[k] for k in range(5) ] ) )
	aPM[3] = jnp.sum( jnp.array( [ Zeta_a4[k]*eta_xi[k] for k in range(5) ] ) )
	aPM[4] = jnp.sum( jnp.array( [ Zeta_a5[k]*eta_xi[k] for k in range(5) ] ) )
	aPM[5] = jnp.sum( jnp.array( [ Zeta_a6[k]*eta_xi[k] for k in range(5) ] ) )
	
	#=======================================
	# Ringdown (RD)
	#=======================================

	af = xi + s4*xi2*eta + s5*xi*eta2 + tt0*xi*eta + 2*jnp.sqrt(3)*eta + t2*eta2 + t3*eta3
	f_RD = func_fRD(af)
	f1 = 0.1*f_RD
	f2 = f_RD

	beta2 = diff_PhaseC_PM(eta,aPM,f2)
	beta1 = PhaseC_PM(eta,aPM,f2) - beta2*f2

	#=======================================
	v = piMf**(1./3)
	Phase_SPA = PhaseC_SPA(eta,aIns,v)
	Phase_PM = PhaseC_PM(eta,aPM,f)
	Phase_RD = PhaseC_RD(f,[beta1,beta2])

	d = 0.005
	PhaseIMR = Phase_SPA*wMinus(f,d,f1) + Phase_PM*wPlus(f,d,f1)*wMinus(f,d,f2) + Phase_RD*wPlus(f,d,f2)

	PhaseC = - PhaseIMR

	#-----------------------------
	# Analytic
	#-----------------------------
	f_max = f2
	v_max = (PI*f_max)**(1./3)

	Phase_SPA_max = PhaseC_SPA(eta,aIns,v_max)
	Phase_PM_max  = PhaseC_PM(eta,aPM,f_max)
	Phase_RD_max  = PhaseC_RD(f_max,[beta1,beta2])

	diff_Phase_SPA_max = diff_PhaseC_SPA(eta,aIns,v_max)
	diff_Phase_PM_max  = diff_PhaseC_PM(eta,aPM,f_max)
	diff_Phase_RD_max  = diff_PhaseC_RD(f_max,[beta1,beta2])

	dwM1 = diff_wMinus(f_max,d,f1) ; wM1 = wMinus(f_max,d,f1)
	dwP1 = diff_wPlus(f_max,d,f1)  ; wP1 = wPlus(f_max,d,f1)
	dwM2 = diff_wMinus(f_max,d,f2) ; wM2 = wMinus(f_max,d,f2)
	dwP2 = diff_wPlus(f_max,d,f2)  ; wP2 = wPlus(f_max,d,f2)

	DiffPhase = diff_Phase_SPA_max*wM1 + Phase_SPA_max*dwM1 \
			  + diff_Phase_RD_max*wP2 + Phase_RD_max*dwP2 \
			  + diff_Phase_PM_max*wP1*wM2 + Phase_PM_max*(dwP1*wM2 + wP1*dwM2)

	PhiC_corr = PhaseC + DiffPhase*f

	return PhiC_corr - PI

def Amplitude_IMRPhenomC(dL,Mc,eta,chi1,chi2,freq):

	m1 = 0.5*(Mc/eta**(3./5)) * (1 + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1 - jnp.sqrt(1-4*eta))

	M = m1+m2
	#eta = m1*m2/M**2
	eta2 = eta*eta
	eta3 = eta2*eta

	GM = M*GMc3
	Mf = freq*GM # dimensionless frequency
	piMf = PI*Mf # := Omega
	V = piMf**(1./3)

	xi = (m1*chi1 + m2*chi2)/M
	xi2 = xi*xi
	xi3 = xi2*xi
	chi_s = (chi1+chi2) # Symmetric Spin Combination
	chi_a = (chi1-chi2) # Anti-Symmetric Spin Combination
	chi_p = chi1*chi2

	Ax = [Mf*(0. + 0.j) for i in range(7)]
	Ax[0] += 1.
	Ax[1] += 0.
	Ax[2] += (-107. + 55.*eta)/42.;
	Ax[3] += 2.*PI - 4.*xi/3. + 2.*eta*chi_s/3.;
	Ax[4] += -2.173/1.512 - eta*(10.69/2.16 - 2.*chi_p) + 2.047*eta2/1.512;
	Ax[5] += -10.7*PI/2.1 + eta*(3.4*PI/2.1) -24.j*eta;
	Ax[6] += 270.27409/6.46800 - 8.56*gE/1.05+ 4.28j*PI/1.05 + 2.*pi2/3. +\
	  eta*(4.1*pi2/9.6 - 27.8185/3.3264) -\
	  20.261*eta2/2.772 + 11.4635*eta3/9.9792 -\
	  4.28*jnp.log(16.*V*V)/1.05;

	Ay = [Mf*0. for i in range(8)]
	Ay[0] += 1.
	Ay[1] += 0.
	Ay[2] += -7.43/3.36 - 11.*eta/4.;
	Ay[3] += 4.*PI - 11.3*xi/1.2 + 19.*eta*chi_s/6.;
	Ay[4] += 3.4103/1.8144 + 5*xi2 + eta*(13.661/2.016 - chi_p/8.) + 5.9*eta2/1.8;
	Ay[5] += -PI*(41.59/6.72 + 189.*eta/8.) - xi*(31.571/1.008 - 116.5*eta/2.4) +\
	      chi_s*(21.863*eta/1.008 - 79.*eta2/6.) - 3*xi*xi2/4. +\
	      9.*eta*xi*chi_p/4.;
	Ay[6] += 164.47322263/1.39708800 - 17.12*gE/1.05 +\
	      16.*pi2/3 - 8.56*jnp.log(16.*V*V)/1.05 +\
	      eta*(45.1*pi2/4.8 - 561.98689/2.17728) +\
	      5.41*eta2/8.96 - 5.605*eta*eta2/2.592 - 80.*PI*xi/3. +\
	      eta*chi_s*(20.*PI/3. - 113.5*xi/3.6) +\
	      xi2*(64.153/1.008 - 45.7*eta/3.6) -\
	      chi_p*(7.87*eta/1.44 - 30.37*eta2/1.44);
	 
	Ay[7] += -PI*(4.415/4.032 - 358.675*eta/6.048 - 91.495*eta2/1.512) -\
	      xi*(252.9407/2.7216 - 845.827*eta/6.048 + 415.51*eta2/8.64) +\
	      chi_s*(158.0239*eta/5.4432 - 451.597*eta2/6.048 + 20.45*eta2*eta/4.32 + 107.*eta*xi2/6. - 5.*eta2*chi_p/24.) +\
	      12.*PI*xi2 - xi2*xi*(150.5/2.4 + eta/8.) +\
	      xi*chi_p*(10.1*eta/2.4 + 3.*eta2/8.);

	SumA = Mf*(0. + 0.j)
	SumX = Mf*0.
	for k in range(7): SumA += Ax[k]*(V**k)
	for k in range(8): SumX += Ay[k]*(V**k)

	R = SumA/jnp.sqrt(jnp.abs(SumX))
	R = jnp.sqrt(R*jnp.conj(R)).real # following line 454 of LALSimIMRPhenomC_internals.c

	#----------------------------------------------------------
	# Computing Coefficients (gamma1, delta1, delta2)
	#----------------------------------------------------------

	eta_xi = [xi,xi2,eta*xi,eta,eta2]

	gamma1 = jnp.sum( jnp.array([ Zeta_g1[k]*eta_xi[k] for k in range(5) ]) )
	delta1 = jnp.sum( jnp.array([ Zeta_d1[k]*eta_xi[k] for k in range(5) ]) )
	delta2 = jnp.sum( jnp.array([ Zeta_d2[k]*eta_xi[k] for k in range(5) ]) )

	# Final Spin: obtained from fit shown in arXiv:0710.3345
	af = xi + s4*xi2*eta + s5*xi*eta2 + tt0*xi*eta + 2*jnp.sqrt(3)*eta + t2*eta2 + t3*eta3
	f_RD = func_fRD(af) # dimensionless
	
	f0 = 0.98*f_RD
	sigma = delta2 * f_RD / func_Q(af)

	d = 0.015
	wM = wMinus(Mf,d,f0)
	wP = wPlus(Mf,d,f0)

	K = jnp.sqrt(1.5/eta)*PI**(1./6)
	A_left = ( R + K*Mf**2*gamma1 )
	A_right = K*delta1*Lorentzian(Mf,f_RD,sigma)
	Coeff = A0 * jnp.sqrt(eta)*M**(5./6)*freq**(-7./6) / dL
	Amp = Coeff * (A_left * wM + A_right * wP)

	return Amp

@jit	
def hphx_IMRPhenomC(dL,iota,phi0,t0,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq): # 12 GwPrms + freq
	chi1, chi2 = sz1, sz2
	Amp   = Amplitude_IMRPhenomC(dL,Mc,eta,chi1,chi2,freq)  # 5 prms + freq
	Phase = Phase_IMRPhenomC(t0,phi0,Mc,eta,chi1,chi2,freq)   # 6 prms + freq

	u = jnp.cos(iota)
	
	gp = 0.5*(1+u*u)
	gx = -1.j*u

	Exp0 = jnp.exp(-1.j*(2*PI*freq*t0 - phi0))

	hp = Amp*gp*jnp.exp(1.j*Phase) * Exp0
	hx = Amp*gx*jnp.exp(1.j*Phase) * Exp0

	return hp, hx