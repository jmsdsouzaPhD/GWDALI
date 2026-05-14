#=========================================================
# IMRPhenomA: IMR waveform for spinless binaries
#=========================================================
# Check LALSimIMRPhenom.c (line 302) 
# [https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_i_m_r_phenom_8c_source.html#l00302]

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

@jit
def Polynomial(a,b,c,eta,M):
	eta2 = eta*eta
	return (a*eta2 + b*eta + c)/(PI*M)

@jit
def Lorentzian(freq,f0,sig):
	sig2 = sig*sig
	df = (freq-f0)
	return 0.5*sig/PI/(df*df + 0.25*sig2)

Psi = []
Psi.append([-1.5829e-1 , 8.7016e-2 , -3.3382e-2 ])
Psi.append([0., 0., 0. ])
Psi.append([3.2967e1 , -1.9000e1 , 2.1345e0])
Psi.append([-3.0849e2 , 1.8211e2 , -2.1727e1])
Psi.append([1.1525e3 , -7.1477e2 , 9.9692e1 ])
Psi.append([0., 0., 0.])
Psi.append([1.2057e3 , -8.4233e2 , 1.8046e2  ])
Psi.append([0. , 0. , 0. ])

@jit
def Phase_IMRPhenomA(t0,phi0,Mc,eta,freq):
	m1 = 0.5*(Mc/eta**(3./5)) * (1. + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1. - jnp.sqrt(1-4*eta))

	M = m1+m2
	M2 = M*M
	eta2 = eta*eta

	GM = M*GMc3
	Mf = GM*freq
	piMf = PI*Mf

	Summation = 0
	for k in range(8):
		xk, yk, zk = Psi[k]
		n = (k-5.)/3
		psi_k = xk*eta2 + yk*eta + zk
		Summation += psi_k*piMf**n

	Phase = Summation/eta - PI 
	return -Phase

A0 = 7.806521525937888e-23

@jit	
def Amplitude_IMRPhenomA(dL,Mc,eta,freq):
	
	m1 = 0.5*(Mc/eta**(3./5)) * (1 + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1 - jnp.sqrt(1-4*eta))

	M = m1+m2
	M2 = M*M
	eta2 = eta*eta

	GM = M*GMc3 

	h0 = A0 * jnp.sqrt(eta)*M**(5./6)/dL

	f_merg  = Polynomial(6.6389e-1,-1.0321e-1,1.0979e-1,eta,GM)
	f_ring  = Polynomial(1.3278,-2.0642e-1,2.1957e-1,eta,GM)
	f_cut   = Polynomial(1.7086,-2.6592e-1,2.8236e-1,eta,GM)
	sig		= Polynomial(1.1383,-1.7700e-1,4.6834e-2,eta,GM)

	w = 0.5*PI*sig*(f_ring/f_merg)**(-2./3)

	ff = [f_merg, f_ring, f_cut, sig]

	Coeff = h0*f_merg**(-7./6)

	A1 = (freq/f_merg)**(-7./6)
	A2 = (freq/f_merg)**(-2./3)
	A3 = w*Lorentzian(freq,f_ring,sig)

	Amp = Coeff*( A1*(freq<f_merg) + (freq>=f_merg)*A2*(freq<f_ring) + (freq>=f_ring)*A3)#*(freq<f_cut) )

	return Amp

@jit	
def hphx_IMRPhenomA(dL,iota,phi0,t0,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq): # 12 GwPrms + freq
	Amp   = Amplitude_IMRPhenomA(dL,Mc,eta,freq)
	Phase = Phase_IMRPhenomA(t0,phi0,Mc,eta,freq)

	u = jnp.cos(iota)
	
	gp = 0.5*(1+u*u)
	gx = -1.j*u

	Exp0 = jnp.exp(-1.j*(2*PI*freq*t0 - phi0))
	hp = Amp*gp*jnp.exp(1.j*Phase) * Exp0
	hx = Amp*gx*jnp.exp(1.j*Phase) * Exp0

	return hp, hx