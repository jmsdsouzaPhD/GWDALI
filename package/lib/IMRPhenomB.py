#=========================================================
# IMRPhenomB: [arXiv:0909.2867v3]
#=========================================================
# Check LALSimIMRPhenom.c (line 302) 
# [https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_i_m_r_phenom_8c_source.html#l00302]

from jax import jit, vmap
import jax.numpy as jnp
import numpy as np

PI = 3.141592653589793238462643383279502884
gE = 0.577215664901532860606512090082402431
Msun = 1.988409870698050731911960804878414216e30
pc = 3.085677581491367278913937957796471611e16
G = 6.67430e-11
c = 299792458.
Mpc = pc*1.e6

GMc3 = 4.925490947641267e-06 # GM/c^3
A0 = 7.806521525937888e-23

pi2 = PI*PI

#===========================================================
# Phase
#===========================================================
@jit
def Phase_IMRPhenomB(t0,phi0,Mc,eta,chi1,chi2,freq):

	m1 = 0.5*(Mc/eta**(3./5)) * (1 + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1 - jnp.sqrt(1-4*eta))

	M = m1+m2
	M2 = M*M
	eta = m1*m2/M2
	eta2 = eta*eta
	eta3 = eta2*eta

	xi = (m1*chi1 + m2*chi2)/M
	xi2 = xi*xi

	GM = PI*M*GMc3
	piMf = GM*freq

	psi0 = [0,0,0,0,0,0,0]
	psi1 = [0,0,0,0,0,0,0]
	psi2 = [  4.914021164021164,		 			-9.2091e2 , 4.9213e2 , 1.3503e2, 6.7419e3, -1.0534e3, -1.3397e4  ]
	psi3 = [-16*PI + 113.*xi/3,					1.7022e4, -9.5659e3, -2.1821e3, -1.2137e5, 2.0752e4, 2.3859e5	]
	psi4 = [30.10315295099521 - 405.*xi2/8, 			-1.2544e5, 7.5066e4, 1.3382e4, 8.7354e5, -1.6573e5, -1.6936e6]
	psi5 = [0,0,0,0,0,0,0]
	psi6 = [0,							-8.8977e5, 6.3102e5, 5.0676e4, 5.9808e6, -1.4148e6, -1.1280e7 ]
	psi7 = [0,							8.6960e5, -6.7098e5, -3.0082e4, -5.8379e6, 1.5145e6, 1.0891e7]
	psi8 = [0,							-3.6600e5, 3.0670e5, 6.3176e2, 2.4265e6, -7.2180e5, -4.5524e6]
	PsiC = [psi0,psi1,psi2,psi3,psi4,psi5,psi6,psi7,psi8]

	etaX = [1., eta, eta*xi, eta*xi2, eta2, eta2*xi, eta3]

	Summation = 0
	for k in range(2,9):
		PsiK = sum( [ PsiC[k][i]*etaX[i] for i in range(7) ] )
		Summation += PsiK*piMf**(k/3)

	Phase = 3./(128*eta*piMf**(5./3))* ( 1. + Summation) - PI
	return -Phase
#---------------------------------------------------------
# Some functions from arXiv:1005.3306v3
#---------------------------------------------------------
@jit
def Lorentzian(f,f0,sig,GM):
	sig2 = sig*sig
	df = (f-f0)*GM
	df2 = df*df
	return sig2/( df2 + sig2/4 )

@jit
def func_Q(a):
	q1, q2, q3 = 0.7, 1.4187, -0.4990
	return q1 + q2*(1.-a)**q3

@jit
def func_fRD(a,GM):
	k1, k2, k3 = 1.5251, -1.1568, 0.1292
	return (k1 + k2*(1.-a)**k3)/(2*PI*GM)

@jit
def get_fprime(freq,f1,GM):
	u1 = (PI*GM*freq)**(1./3)
	u2 = u1*u1
	u3 = u2*u1
	f_prime = freq/f1
	return u1, u2, u3, f_prime

@jit
def IMRB1(freq,f1,alpha2,alpha3,GM):
	u1, u2, u3, f_prime = get_fprime(freq,f1,GM)
	imr = f_prime**(-7./6)*(1 + alpha2*u2 + alpha3*u3)
	return imr

@jit
def IMRB2(freq,f1,eps1,eps2,GM):
	u1, u2, u3, f_prime = get_fprime(freq,f1,GM)
	imr = f_prime**(-2./3)*(1 + eps1*u1 + eps2*u2)
	return imr

@jit
def IMRB3(freq,f1,f2,sigma,GM):
		imr = Lorentzian(freq,f2,sigma*GM,GM)
		return imr

@jit
def Amplitude_IMRPhenomB(dL,Mc,eta,chi1,chi2,freq):

	m1 = 0.5*(Mc/eta**(3./5)) * (1 + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1 - jnp.sqrt(1-4*eta))

	M = m1+m2
	eta2 = eta*eta
	eta3 = eta2*eta
	GM = M*GMc3

	C = A0 * Mc**(5./6) / dL 

	a = chi = (m1*chi1 + m2*chi2)/M
	a2 = a*a
	
	dXi = 1-chi
	mu01 = 1.-4.4547*dXi**0.217+3.521*dXi**0.26
	mu02 = 0.5*( 1. - 0.63*dXi**0.3 )
	mu0S = 0.25*( 1. - 0.63*dXi**0.3 )*dXi**0.45
	mu03 = 3.2361e-1 + 4.8935e-2*chi + 1.3463e-2*chi*chi

	row1 = [mu01, 6.4365e-1, 8.2696e-1, -2.7063e-1, -5.8218e-2, -3.9346, -7.0916 ]
	row2 = [mu02, 1.4690e-1,-1.2281e-1,-2.6091e-2,-2.4900e-2, 1.7013e-1, 2.3252 ]
	rowS = [mu0S, -4.0979e-1,-3.5226e-2,1.0082e-1,1.8286,-2.0169e-2,-2.8698]
	row3 = [mu03, -1.3313e-1,-8.1719e-2,1.4512e-1,-2.7140e-1,1.2788e-1,4.9220]

	Coeff = 64625.006840429545 / M
	f1    = Coeff*( row1[0] + row1[1]*eta + row1[2]*eta*chi + row1[3]*eta*chi*chi + row1[4]*eta2 + row1[5]*eta2*chi + row1[6]*eta3 )
	f2    = Coeff*( row2[0] + row2[1]*eta + row2[2]*eta*chi + row2[3]*eta*chi*chi + row2[4]*eta2 + row2[5]*eta2*chi + row2[6]*eta3 )
	sigma = Coeff*( rowS[0] + rowS[1]*eta + rowS[2]*eta*chi + rowS[3]*eta*chi*chi + rowS[4]*eta2 + rowS[5]*eta2*chi + rowS[6]*eta3 )
	f3    = Coeff*( row3[0] + row3[1]*eta + row3[2]*eta*chi + row3[3]*eta*chi*chi + row3[4]*eta2 + row3[5]*eta2*chi + row3[6]*eta3 )

	alpha2 = -323./224 + 451*eta/168
	alpha3 = (27./8-11*eta/6)*chi
	eps1 = 1.4547*chi-1.8897
	eps2 = -1.8153*chi + 1.6557

	Const = C*f1**(-7/6)

	AmpB_1  = IMRB1(freq,f1,alpha2,alpha3,GM)
	imr1_f1 = IMRB1(f1,f1,alpha2,alpha3,GM)
	imr2_f1 = IMRB2(f1,f1,eps1,eps2,GM)
	wm = imr1_f1/imr2_f1
	imr2_f2 =  wm*IMRB2(f2,f1,eps1,eps2,GM)
	imr3_f2 = IMRB3(f2,f1,f2,sigma,GM)
	wr = imr2_f2/imr3_f2

	AmpB_2 = wm*IMRB2(freq,f1,eps1,eps2,GM)
	AmpB_3 = wr*IMRB3(freq,f1,f2,sigma,GM)

	Amp = Const * ( AmpB_1*(freq<f1) + (freq>=f1)*AmpB_2*(freq<f2) + (freq>=f2)*AmpB_3)#*(freq<f3) )

	return Amp

@jit
def hphx_IMRPhenomB(dL,iota,phi0,t0,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq): # 12 GwPrms + freq
	chi1, chi2 = sz1, sz2
	Amp   = Amplitude_IMRPhenomB(dL,Mc,eta,chi1,chi2,freq)  # 5 prms + freq
	Phase = Phase_IMRPhenomB(t0,phi0,Mc,eta,chi1,chi2,freq)   # 6 prms + freq

	u = jnp.cos(iota)
	
	gp = 0.5*(1+u*u)
	gx = -1.j*u

	Exp0 = jnp.exp(-1.j*(2*PI*freq*t0 - phi0))
	hp = Amp*gp*jnp.exp(1.j*Phase) * Exp0
	hx = Amp*gx*jnp.exp(1.j*Phase) * Exp0

	return hp, hx