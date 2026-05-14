import jax
import jax.numpy as jnp
from jax import jit, vmap

import numpy as np
import pandas as pd
from pathlib import Path

PI = 3.141592653589793238462643383279502884
gE = 0.577215664901532860606512090082402431
Msun = 1.988409870698050731911960804878414216e30
pc = 3.085677581491367278913937957796471611e16
G = 6.67430e-11
c = 299792458.
Mpc = pc*1.e6
Gpc = pc*1.e9

PI_M_SIXTH = 0.8263074871107581108331125856317241299
GMc3 = G*Msun/c**3 # Equal to LAL_MTSUN_SI
#GMc3 = 4.925490947641266978197229498498379006e-6

A0 = 7.806521525937888e-23
pi2 = PI*PI

path = Path(__file__).parent / '.'
Table_Lambda = pd.read_csv(path / 'table.csv')
data = np.loadtxt(path / 'QNMData.txt')
aa  = data[:,0]
fRg = data[:,1]
fDp = data[:,2]

data_rows = [ Table_Lambda.iloc[i] for i in range(19) ]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# JAX CubicSpline Implementation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@jit
def compute_cubic_spline_coefficients_jax(x, y):
	x = jnp.array(x,dtype=jnp.float64)
	y = jnp.array(y,dtype=jnp.float64)

	a = y.astype(jnp.float64).copy()
	h = jnp.diff(x)
	alpha = jnp.zeros(len(x))

	def compute_alpha(i, alpha):
		return alpha.at[i].set((3/h[i]) * (a[i+1] - a[i]) - (3/h[i-1]) * (a[i] - a[i-1]))

	alpha = jax.lax.fori_loop(1, len(x)-1, compute_alpha, alpha)
	
	l  = jnp.ones(len(x)-1+1 , dtype=jnp.float64)
	mu = jnp.zeros(len(x)-1+1 , dtype=jnp.float64)
	z  = jnp.zeros(len(x)-1+1 , dtype=jnp.float64)
	c  = jnp.zeros(len(x)-1+1 , dtype=jnp.float64)
	b  = jnp.zeros(len(x)-1 , dtype=jnp.float64)
	d  = jnp.zeros(len(x)-1 , dtype=jnp.float64)

	def compute_l_mu_z(i, l_mu_z):
		l, mu, z = l_mu_z
		l = l.at[i].set(2 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1])
		mu = mu.at[i].set(h[i] / l[i])
		z = z.at[i].set((alpha[i] - h[i-1] * z[i-1]) / l[i])
		return (l, mu, z)#, None

	(l, mu, z) = jax.lax.fori_loop(1, len(x)-1, compute_l_mu_z, (l, mu, z))
	

	l = l.at[len(x)-1].set(1)
	z = z.at[len(x)-1].set(0)
	c = c.at[len(x)-1].set(0)

	def compute_c_b_d(c_b_d,j):
		c, b, d = c_b_d
		c = c.at[j].set(z[j] - mu[j] * c[j+1])
		b = b.at[j].set((a[j+1] - a[j]) / h[j] - h[j] * (c[j+1] + 2 * c[j]) / 3)
		d = d.at[j].set((c[j+1] - c[j]) / (3 * h[j]))
		return (c, b, d),None

	(c, b, d), _ = jax.lax.scan(compute_c_b_d, (c, b, d), jnp.arange(len(x)-2, -1, -1))
	return a, b, c, d

@jit
def cubic_spline_evaluate_jax(x_new, x, a, b, c, d):
    i = jnp.searchsorted(x, x_new) - 1
    i = jnp.clip(i, 0, len(x) - 2)
    dx = x_new - x[i]
    return a[i] + b[i] * dx + c[i] * dx**2 + d[i] * dx**3

def JAX_CubicSpline(x_new, x, y):
    a, b, c, d = compute_cubic_spline_coefficients_jax(x, y)
    return cubic_spline_evaluate_jax(x_new, x, a, b, c, d)

func_interp_RD = lambda x_new: jnp.interp(x_new,aa,fRg)
func_interp_damp = lambda x_new: jnp.interp(x_new,aa,fDp)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@jit
def Final_Spin(eta,S,S_hat):
	eta2 = eta*eta
	eta3 = eta2*eta
	eta4 = eta3*eta

	#================================================
	# Final Spin: Eq. 3.6 of arXiv:1508.07250
	#================================================
	S2 = S*S
	S3 = S2*S
	S4 = S3*S

	af = eta*(3.4641016151377544 - 4.399247300629289*eta + \
	      9.397292189321194*eta2 - 13.180949901606242*eta3 + \
	      S*((1.0/eta - 0.0850917821418767 - 5.837029316602263*eta) + \
	      (0.1014665242971878 - 2.0967746996832157*eta)*S + \
	      (-1.3546806617824356 + 4.108962025369336*eta)*S2 + \
	      (-0.8676969352555539 + 2.064046835273906*eta)*S3))

	#================================================
	# Radieted Energy: Eqs. 3.7 and 3.8 of arXiv:1508.07250
	#================================================

	E_rad = eta * (0.055974469826360077 + 0.5809510763115132 * eta - 0.9606726679372312 * eta2 + 3.352411249771192 * eta3) * \
	            (1. + (-0.0030302335878845507 - 2.0066110851351073 * eta + 7.7050567802399215 * eta2) * S_hat) / \
	           (1. + (-0.6714403054720589 - 1.4756929437702908 * eta + 7.304676214885011 * eta2) * S_hat)

	return af, E_rad

#==============================================================#==============================================================
@jit
def PhiD_TF2(Mf,eta,varphi):
	# LALSimIMRPhenomD_internals.c does not include (2*pi*f*t_c - phi_c) at TaylorF2 Phase ("PhiInsAnsatzInt")
	# The term (2*pi*f*t_c - phi_c) enters only in "Phase_IMRPhenomD(m1,m2,chi1,chi2,t0,phi0,freq)" bellow
	piMf = PI*Mf
	v = (piMf)**(1./3)
	v5 = v**5
	ones = Mf*0 + 1
	VarPhi = [x*ones for x in varphi]
	# [5]: f-dependent varphi[5]*=(1+jnp.log(piMf))
	# [6]: f-dependent: varphi[6]-=6848*jnp.log(64*piMf)/63
	VarPhi[5] *= (1+jnp.log(piMf))
	VarPhi[6] -= 6848*jnp.log(64*piMf)/63
	Summation = Mf*0
	for k in range(8):
		Summation += VarPhi[k]*v**k
	return 3./(128*eta*v5)*Summation - PI/4

@jit
def PhiD_Ins(Mf,eta,varphi,Sigma):
	phi_tf2 = PhiD_TF2(Mf,eta,varphi)
	sigma1, sigma2, sigma3, sigma4 = Sigma
	# From Eq (28) of IMRPhenomD paper2
	inv_eta = 1./eta
	return phi_tf2 + inv_eta*(sigma1*Mf + 0.75*sigma2*Mf**(4./3) + 0.6*sigma3*Mf**(5./3) + 0.5*sigma4*Mf*Mf)

@jit
def PhiD_Int(Mf,eta,Beta):
	beta1, beta2, beta3 = Beta
	inv_eta = 1./eta
	return inv_eta * ( beta1*Mf + beta2*jnp.log(Mf) - beta3/(3*Mf**3) )

@jit
def PhiD_MR(Mf,eta,Alpha,f_RD,f_damp):
	alpha1, alpha2, alpha3, alpha4, alpha5 = Alpha
	df = Mf-alpha5*f_RD
	inv_eta = 1./eta
	return inv_eta * ( alpha1*Mf -alpha2/Mf + (4./3)*alpha3*Mf**(3./4) + alpha4*jnp.arctan( df/f_damp ) )

#==============================================================#==============================================================
@jit
def Diff_PhiD_TF2(Mf,eta,varphi):
	piMf = PI*Mf
	v = (piMf)**(1./3)
	v5 = v**5
	ones = Mf*0+1
	VarPhi = [x*ones for x in varphi]
	# [5]: f-dependent varphi[5]*=(1+jnp.log(piMf))
	# [6]: f-dependent: varphi[6]-=6848*jnp.log(64*piMf)/63
	VarPhi[5] *= (1+jnp.log(piMf))
	VarPhi[6] -= 6848*jnp.log(64*piMf)/63

	diff_VarPhi = [Mf*0 for i in range(len(varphi))]
	diff_VarPhi[5] = varphi[5]*ones / Mf
	diff_VarPhi[6] = -(6848./63) / Mf
	Summation = Mf*0
	for k in range(8):
		Summation += ( (k-5)*VarPhi[k] / 3/Mf + diff_VarPhi[k] ) * v**k
	return 3*Summation/(128*eta*v5)

@jit
def Diff_PhiD_Ins(Mf,eta,varphi,Sigma):
	diff_phi_tf2 = Diff_PhiD_TF2(Mf,eta,varphi)
	sigma1, sigma2, sigma3, sigma4 = Sigma
	return diff_phi_tf2 + (sigma1 + sigma2*Mf**(1./3) + sigma3*Mf**(2./3)+sigma4*Mf)/eta

@jit
def Diff_PhiD_Int(Mf,eta,Beta):
	beta1, beta2, beta3 = Beta
	inv_eta = 1./eta
	return inv_eta*( beta1 + beta2/Mf + beta3/Mf**4 )

@jit
def Diff_PhiD_MR(Mf,eta,Alpha,f_RD,f_damp):
	alpha1, alpha2, alpha3, alpha4, alpha5 = Alpha
	df = Mf-alpha5*f_RD ; df2 = df*df
	f_damp2 = f_damp*f_damp
	inv_eta = 1./eta
	return inv_eta*( alpha1 + alpha2/Mf**2 + alpha3*Mf**(-1/4) + alpha4*f_damp/(f_damp2+df2) )
#==============================================================#==============================================================

#@jit
def Lambda(i,DataFrame,eta,dX):
	dX2 = dX*dX
	dX3 = dX2*dX
	eta2 = eta*eta
	row = DataFrame.iloc[i]
	l = list(row)[1:]
	Lambda_i = l[0] + l[1]*eta + \
			   dX *(l[2]+l[3]*eta+l[4]*eta2) + \
			   dX2*(l[5]+l[6]*eta+l[7]*eta2) + \
			   dX3*(l[8]+l[9]*eta+l[10]*eta2)
	return Lambda_i

@jit
def Phase_Connections(f1,f2,eta,varphi,Alpha,Beta,Sigma,f_RD,f_damp):
	# C0, C1 : Coeffitients for Continuity between Phase_Ins and Phase_Int at f=f1:
	C1_int = Diff_PhiD_Ins(f1,eta,varphi,Sigma) - Diff_PhiD_Int(f1,eta,Beta)
	C0_int = PhiD_Ins(f1,eta,varphi,Sigma) - ( PhiD_Int(f1,eta,Beta) + C1_int*f1 )
	
	# C0 and C1 : Coeffitients for Continuity between Phase_Int and Phase MR at f=f2:
	C1_mr = (Diff_PhiD_Int(f2,eta,Beta)+C1_int) - Diff_PhiD_MR(f2,eta,Alpha,f_RD,f_damp)
	C0_mr = (PhiD_Int(f2,eta,Beta)+C0_int+C1_int*f2) - ( PhiD_MR(f2,eta,Alpha,f_RD,f_damp) + C1_mr*f2 )

	return C0_int, C1_int, C0_mr,C1_mr

@jit
def VarPhi_InsPhase(m1,m2,chi1,chi2):
	M = m1+m2 ; M2 = M*M

	m1M = m1/M ; m2M = m2/M

	eta = m1*m2/M2
	eta2 = eta*eta
	eta3 = eta2*eta

	deltaM = m1M - m2M
	chi_s = (chi1+chi2)/2
	chi_a = (chi1-chi2)/2

	xi_eff = chi1*m1M + chi2*m2M
	xi2 = xi_eff*xi_eff
	xi3 = xi2*xi_eff

	varphi = [0]*8

	varphi[0] += 1.
	varphi[2] += (3715./756.) + 55.*eta/9
	varphi[3] += -16.*PI
	varphi[4] += 15293365./508032 + 27145*eta/504 + 3085*eta2/72
	# [5]: f-dependent varphi[5]*=(1+np.log(piMf))
	varphi[5] += 38645*PI/756 - 65*PI*eta/9
	# [6]: f-dependent: varphi[6]-=6848*jnp.log(64*piMf)/63
	varphi[6] += 11583231236531./4694215680 - 6848.*gE/21 - 640*pi2/3 + (-15737765635./3048192 \
				+ 2255*pi2/12)*eta + 76055*eta2/1728 - 127825*eta3/1296
	varphi[7] += 77096675*PI/254016 + 378515*PI*eta/1512 - 74045*PI*eta2/756
	#======================================================================
	# Spinning Corrections to varphi:
	varphi[3] += 113.*deltaM*chi_a/3 + (113./3 - 76.*eta/3)*chi_s
	varphi[4] += (-405./8 + 200*eta)*chi_a**2 - 405*deltaM*chi_a*chi_s/4 + (-405./8 + 5*eta/2)*chi_s**2
	varphi[5] += deltaM*(-732985./2268 - 140*eta/9)*chi_a + (-732985./2268 + 24260*eta/81 + 340*eta2/9)*chi_s
	varphi[6] += 2270*PI*deltaM*chi_a/3 + (2270*PI/3 - 520*PI*eta)*chi_s
	varphi[7] += deltaM*(-25150083775./3048192 + 26804935.*eta/6048-1985*eta2/48)*chi_a + (-25150083775./3048192 + 10566655595.*eta/762048 - 1042165.*eta2/3024 + 5345*eta3/36)*chi_s
	
	return varphi

@jit
def Phase_IMR_D(m1,m2,chi1,chi2,Coeffs,Connections,varphi,f1,f2,f_RD,f_damp,Mf):
	# Coeffs = [g2,g3,s1,s2,s3,s4,b1,b2,b3,a1,a2,a3,a4]
	Sigma = Coeffs[2:6]
	Beta  = Coeffs[6:9]
	Alpha = Coeffs[9:]

	C0_int, C1_int, C0_mr,C1_mr = Connections

	eta = m1*m2/(m1+m2)**2

	Phase_Ins = PhiD_Ins(Mf,eta,varphi,Sigma)
	Phase_Int = PhiD_Int(Mf,eta,Beta) + C0_int + C1_int*Mf
	Phase_MR  = PhiD_MR(Mf,eta,Alpha,f_RD,f_damp) + C0_mr + C1_mr*Mf

	Phase_IMR = Phase_Ins*(Mf<f1) + (Mf>=f1)*Phase_Int*(Mf<f2) + (Mf>=f2)*Phase_MR
	
	return Phase_IMR

@jit
def Phase_IMRPhenomD(*args_phase):
	Mc, eta, chi1, chi2, t0, phi0, f_RD, f_damp, freq = args_phase

	m1 = 0.5*(Mc/eta**(3./5)) * (1 + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1 - jnp.sqrt(1-4*eta))

	M = m1 + m2 ; Ms = M*M
	m1s = m1*m1 
	m2s = m2*m2
	eta = m1*m2/Ms
	GM = M*GMc3
	Mf = GM*freq
	fRef = 1.0
	Mf_ref = fRef*GM

	f1 = 0.018
	f2 = 0.5*f_RD

	#-------------------------------------------------------------
	# Computing Alpha,Beta,Gamma,Sigma Coefficients:
	#-------------------------------------------------------------
	xi_eff = (chi1*m1 + chi2*m2)/M
	xi_PN = xi_eff - (38./113.)*eta*(chi1+chi2)

	dXi = xi_PN - 1.
	# Coeffs = [gamma2,gamma3,
	#			sigma1,sigma2,sigma3,sigma4,
	#			beta1,beta2,beta3,
	#			alpha1,alpha2,alpha3,alpha4,alpha5]
	Coeffs = [ Lambda(i,  Table_Lambda, eta, dXi) for i in range(5,19)]
	gamma2, gamma3 = Coeffs[:2]
	Sigma = Coeffs[2:6]
	Beta  = Coeffs[6:9]
	Alpha = Coeffs[9:]
	#-------------------------------------------------------------
	varphi = VarPhi_InsPhase(m1,m2,chi1,chi2)

	Connections = Phase_Connections(f1,f2,eta,varphi,Alpha,Beta,Sigma,f_RD,f_damp)

	Phase   = Phase_IMR_D(m1,m2,chi1,chi2,Coeffs,Connections,varphi,f1,f2,f_RD,f_damp,Mf)
	phi_ref = Phase_IMR_D(m1,m2,chi1,chi2,Coeffs,Connections,varphi,f1,f2,f_RD,f_damp,Mf_ref)
	
	#----------------------------------------------------------------------------
	C0_int, C1_int, C0_mr,C1_mr = Connections
	alpha5 = Alpha[4]
	gamma22 = gamma2*gamma2
	f_peak = jnp.abs( f_RD + (f_damp*gamma3/gamma2 )*(jnp.sqrt(1.-gamma22)-1.))
	df     = f_peak-alpha5*f_RD ; df2 = df*df
	f_damp2 = f_damp*f_damp
	t_corr = Diff_PhiD_MR(f_peak,eta,Alpha,f_RD,f_damp)
	#----------------------------------------------------------------------------

	Phase -= (Mf-Mf_ref)*t_corr + phi_ref - (2*PI*freq*t0 - phi0)
	return -Phase

#=====================================================#=====================================================
#------------------------------------------------Amplitude--------------------------------------------------
#=====================================================#=====================================================

@jit
def funcA0(x,eta):
	return jnp.sqrt(eta)*x**(-7./6)

#=====================================================
# Inspiral (Function and Derivative)
#=====================================================

@jit
def Amplitude_Ins(x,eta,Ai,rho):
	A_PN = 0 ; A_rho = 0
	K = funcA0(x,eta)
	for i in range(7): A_PN += Ai[i]*(PI*x)**(i/3)
	for i, rho_i in enumerate(rho): A_rho += rho_i*x**((7.+i)/3) 

	A_Ins = K*(A_PN + A_rho)
	return A_Ins

@jit
def diff_Ins(x,eta,Ai,rho):
	Ains = Amplitude_Ins(x,eta,Ai,rho)
	SumA   = 0 
	SumRho = 0
	K = funcA0(x,eta)
	for k in range(7):
		SumA += k*Ai[k]*(PI*x)**(k/3)
	for k, rho_i in enumerate(rho):
		SumRho += rho_i*(7.+k)*x**((7.+k)/3)
	return ( -7*Ains/2 + K*(SumA+SumRho) )/(3*x)

#=====================================================
# Merger-Ringdown (Function and Derivative)
#=====================================================
@jit
def Amplitude_MR(x,eta,f_qnm,gamma):
	f_RD, f_damp = f_qnm
	gamma1, gamma2, gamma3 = gamma
	df = (x-f_RD)
	df2 = df*df
	arg = -gamma2*df/(gamma3*f_damp)
	K = funcA0(x,eta)
	A_MR = K*gamma1*gamma3*f_damp/(df2+(gamma3*f_damp)**2) * jnp.exp(arg)
	return A_MR

@jit
def diff_MR(x,eta,f_qnm,gamma):
	f_RD, f_damp = f_qnm
	gamma1, gamma2, gamma3 = gamma
	df = x-f_RD
	df2 = df*df
	D = df2+(gamma3*f_damp)**2
	Amr = Amplitude_MR(x,eta,f_qnm,gamma)
	return -Amr*( (7*gamma3*f_damp + 6*gamma2*x)*D +12*gamma3*f_damp*df*x)/(6*gamma3*f_damp*x*D)

#=====================================================
# Intermediate (Function and Derivative)
#=====================================================
@jit
def Amplitude_Int(x,eta,delta):
	del0,del1,del2,del3,del4 = delta
	K = funcA0(x,eta)
	return K*( del0 + del1*x + del2*x**2 + del3*x**3 + del4*x**4 )

@jit
def diff_Int(x,eta,delta):
	del0,del1,del2,del3,del4 = delta
	K = funcA0(x,eta)
	return (-7*del0 - x*del1 + 5*del2*x**2 + 11*del3*x**3 + 17*del4*x**4 ) *K/(6*x)

@jit
def Amplitude_IMRPhenomD(*args_amp):
	dL, iota, Mc, eta, chi1, chi2, f_RD, f_damp, freq = args_amp

	m1 = 0.5*(Mc/eta**(3./5)) * (1 + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1 - jnp.sqrt(1-4*eta))
	m1s, m2s = m1*m1, m2*m2

	M = m1+m2
	delta = (m1-m2)/M
	eta2 = eta*eta
	eta3 = eta2*eta
	eta4 = eta3*eta

	GM = M*GMc3
	f = freq*GM

	chi_eff = (chi1*m1 + chi2*m2)/M
	chi_PN  = chi_eff - 38*eta*(chi1+chi2)/113
	chi_s = (chi1+chi2)/2
	chi_a = (chi1-chi2)/2
	chi_eff2 = chi_eff**2
	chi_eff3 = chi_eff**3
	chi_eff4 = chi_eff**4

	f_qnm = [f_RD, f_damp]

	#================================================
	# Coefficients
	#================================================

	# From Eqs. (B14) - (B20) of [PhysRevD.93.044007]

	Ai = [0]*7
	Ai[0] = 1.
	Ai[1] = 0.
	Ai[2] = -(323./224) + (451./168)*eta
	Ai[3] = 27*delta*chi_a/8 + (27./8 - 11*eta/6)*chi_s
	Ai[4] = -27312085./8128512 - 1975055*eta/338688 + 105271*eta2/24192 +\
	 (-81./32 + 8*eta)*chi_a**2 - 81.*delta*chi_a*chi_s/16 + (-81./32 + 17*eta/8)*chi_s**2
	Ai[5] = -85*PI/64 + 85*PI*eta/16 + delta*(285197./16128 - 1579*eta/4032)*chi_a \
	+ (285197./16128 - 15317*eta/672 - 2227*eta2/1008)*chi_s
	Ai[6] = -177520268561./8583708672 + (545384828789./5007163392 - 205*pi2/48)*eta \
	- 3248849057.*eta2/178827264 + 34473079*eta3/6386688 \
	+ (1614569./64512 - 1873643*eta/16128 + 2167*eta2/42)*chi_a**2 \
	+ (31*PI/12 - 7*PI*eta/3)*chi_s + (1614569./64512 - 61391*eta/1344 + 57451*eta2/4032)*chi_s**2 \
	+ delta*chi_a*(31*PI/12 + (1614569./32256 - 165961*eta/2688)*chi_s)

	dXi = chi_PN - 1
	rho1   = Lambda(0, Table_Lambda, eta, dXi)
	rho2   = Lambda(1, Table_Lambda, eta, dXi)
	rho3   = Lambda(2, Table_Lambda, eta, dXi)
	u2     = Lambda(3, Table_Lambda, eta, dXi)
	gamma1 = Lambda(4, Table_Lambda, eta, dXi)
	gamma2 = Lambda(5, Table_Lambda, eta, dXi)
	gamma3 = Lambda(6, Table_Lambda, eta, dXi)

	rho = [rho1,rho2,rho3]
	gamma = [gamma1, gamma2, gamma3]

	f1 = 0.014
	f3 = jnp.abs( f_RD + f_damp*gamma3*(jnp.sqrt(1.-gamma2**2)-1.)/gamma2 )
	f2 = (f1+f3)/2

	u1 = Amplitude_Ins(f1,eta,Ai,rho)
	u3 = Amplitude_MR(f3,eta,f_qnm,gamma)
	d1 = diff_Ins(f1,eta,Ai,rho)
	d3 = diff_MR(f3,eta,f_qnm,gamma)

	K1 = funcA0(f1,eta)
	K2 = funcA0(f2,eta)
	K3 = funcA0(f3,eta)
	u2 *= K2

	fpoly = lambda ff,a: jnp.array([a[n]*ff**n for n in range(5)])

	r1 = fpoly(f1,jnp.ones(5))*K1
	r2 = fpoly(f2,jnp.ones(5))*K2
	r3 = fpoly(f3,jnp.ones(5))*K3
	r4 = fpoly(f1,[-7.,-1.,5.,11.,17.])*K1 / (6.*f1)
	r5 = fpoly(f3,[-7.,-1.,5.,11.,17.])*K3 / (6.*f3)

	Matrix = jnp.array([r1,r2,r3,r4,r5])
	delta = jnp.linalg.solve(Matrix,jnp.array([u1,u2,u3,d1,d3]))
	
	A_Ins = Amplitude_Ins(f,eta,Ai,rho)
	A_MR  = Amplitude_MR(f,eta,f_qnm,gamma)
	A_Int = Amplitude_Int(f,eta,delta)

	A_IMR = A_Ins*(f<f1) + (f>=f1)*A_Int*(f<f3) + (f>=f3)*A_MR

	Amp = A0 * A_IMR * M*M*(GMc3**(7./6))/dL

	return Amp

@jit	
def hphx_IMRPhenomD(dL,iota,phi0,t0,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq): # 12 GwPrms + freq
	chi1, chi2 = sz1, sz2
	#================================================
	# Ringdown and Damping Frequencies
	#================================================
	m1 = 0.5*(Mc/eta**(3./5)) * (1 + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1 - jnp.sqrt(1-4*eta))
	M = m1+m2
	
	S = (m1*m1*chi1 + m2*m2*chi2)/M**2
	S_hat = S_hat = (m1*m1*chi1 + m2*m2*chi2)/(m1*m1+m2*m2)

	af, E_rad = Final_Spin(eta,S,S_hat)
	
	f_RD   = func_interp_RD(af) / (1.-E_rad)  
	f_damp = func_interp_damp(af) / (1.-E_rad)
	
	#f_RD   = JAX_CubicSpline(af, aa, fRg) / (1.-E_rad) # Recovery lalsimulation
	#f_damp = JAX_CubicSpline(af, aa, fDp) / (1.-E_rad) # Recovery lalsimulation

	#================================================

	args_amp   = [dL, iota, Mc, eta, chi1, chi2, f_RD, f_damp, freq]
	args_phase = [Mc, eta, chi1, chi2, t0, phi0, f_RD, f_damp, freq]

	Amp   = Amplitude_IMRPhenomD(*args_amp)
	Phase = Phase_IMRPhenomD(*args_phase)

	u = jnp.cos(iota)
	
	gp = 0.5*(1+u*u)
	gx = -1.j*u

	hp = Amp*gp*jnp.exp(1.j*Phase)
	hx = Amp*gx*jnp.exp(1.j*Phase)

	return hp, hx