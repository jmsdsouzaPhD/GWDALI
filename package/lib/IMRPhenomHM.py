import jax
import jax.numpy as jnp
from jax import jit, vmap, debug

#jax.config.update("jax_disable_jit",True) # Desabilitar jax.jit

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

idxs_HM = [[2,2],[2,1],[3,3],[3,2],[4,4],[4,3]]

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
def PhiD_TF2(Mf,eta,varphi):
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

	Summation = jnp.array(Summation,dtype=jnp.float64)
	VarPhi = jnp.array(VarPhi,dtype=jnp.float64)
	v = jnp.array(v,dtype=jnp.float64)
	def add_varphi(k,Summation):
		return Summation + VarPhi[k]*v**k
	Summation = jax.lax.fori_loop(0, 8, add_varphi, Summation)
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
def PhiD_MR(Mf,eta,Alpha,f_RD,f_damp, rho_lm, tau_lm):
	alpha1, alpha2, alpha3, alpha4, alpha5 = Alpha
	df = Mf-alpha5*f_RD
	fdp = rho_lm*tau_lm*f_damp
	alph4 = alpha4*rho_lm
	inv_eta = 1./eta
	return inv_eta * ( alpha1*Mf -alpha2/Mf + (4./3)*alpha3*Mf**(3./4) + alph4*jnp.arctan( df/fdp ) )

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
	
	Summation = jnp.array(Summation,dtype=jnp.float64)
	VarPhi = jnp.array(VarPhi,dtype=jnp.float64)
	diff_VarPhi = jnp.array(diff_VarPhi,dtype=jnp.float64)
	v = jnp.array(v,dtype=jnp.float64)
	Mf = jnp.array(Mf,dtype=jnp.float64)
	def add_diff_vaprhi(k,Summation):
		return Summation + ( (k-5)*VarPhi[k] / 3/Mf + diff_VarPhi[k] ) * v**k
	Summation = jax.lax.fori_loop(0,8,add_diff_vaprhi,Summation)
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
def Diff_PhiD_MR(Mf,eta,Alpha,f_RD,f_damp, rho_lm, tau_lm):
	alpha1, alpha2, alpha3, alpha4, alpha5 = Alpha
	df = Mf-alpha5*f_RD ; df2 = df*df
	fdp = f_damp*(rho_lm*tau_lm)
	alph4 = alpha4*rho_lm
	fdp2 = fdp*fdp
	inv_eta = 1./eta
	return inv_eta*( alpha1 + alpha2/Mf**2 + alpha3*Mf**(-1./4) + alph4*fdp/(fdp2+df2) )
#==============================================================#==============================================================

Lambda_Matrix = jnp.array(Table_Lambda.iloc[0:,1:].to_numpy(), dtype=jnp.float64)

@jit
def Lambda(i,eta,dX):
	dX2 = dX*dX
	dX3 = dX2*dX
	eta2 = eta*eta
	row = jnp.array(Lambda_Matrix[i], dtype=jnp.float64)
	Lambda_i = row[0] + row[1]*eta + \
			   dX  * (row[2] + row[3]*eta + row[4]*eta2) + \
			   dX2 * (row[5] + row[6]*eta + row[7]*eta2) + \
			   dX3 * (row[8] + row[9]*eta + row[10]*eta2)
	return Lambda_i

@jit
def Phase_Connections(f1,f2,eta,varphi,Alpha,Beta,Sigma,f_RD,f_damp,rho_lm,tau_lm):
	# C0, C1 : Coeffitients for Continuity between Phase_Ins and Phase_Int at f=f1:
	C1_int = Diff_PhiD_Ins(f1,eta,varphi,Sigma) - Diff_PhiD_Int(f1,eta,Beta)
	C0_int = PhiD_Ins(f1,eta,varphi,Sigma) - ( PhiD_Int(f1,eta,Beta) + C1_int*f1 )
	
	# C0 and C1 : Coeffitients for Continuity between Phase_Int and Phase MR at f=f2:
	C1_mr = (Diff_PhiD_Int(f2,eta,Beta)+C1_int) - Diff_PhiD_MR(f2,eta,Alpha,f_RD,f_damp,rho_lm,tau_lm)
	C0_mr = (PhiD_Int(f2,eta,Beta)+C0_int+C1_int*f2) - ( PhiD_MR(f2,eta,Alpha,f_RD,f_damp,rho_lm,tau_lm) + C1_mr*f2 )

	return C0_int, C1_int, C0_mr,C1_mr

@jit
def VarPhi_InsPhase(eta,deltaM,chi_a,chi_s):
	eta2 = eta*eta
	eta3 = eta2*eta

	varphi = [0]*8
	
	varphi[0] += 1.
	varphi[2] += (3715./756.) + 55.*eta/9
	varphi[3] += -16.*PI
	varphi[4] += 15293365./508032 + 27145*eta/504 + 3085*eta2/72
	# [5]: f-dependent varphi[5]*=(1+jnp.log(piMf))
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
def Phase_IMR_D(eta,Coeffs,Connections,varphi,f1,f2,f_RD,f_damp,rho_lm,tau_lm,Mf):
	Sigma = Coeffs[7:11]
	Beta  = Coeffs[11:14]
	Alpha = Coeffs[14:19]

	C0_int, C1_int, C0_mr,C1_mr = Connections

	Phase_Ins = PhiD_Ins(Mf,eta,varphi,Sigma)
	Phase_Int = PhiD_Int(Mf,eta,Beta) + C0_int + C1_int*Mf
	Phase_MR  = PhiD_MR(Mf,eta,Alpha,f_RD,f_damp,rho_lm,tau_lm) + C0_mr + C1_mr*Mf

	Phase_IMR = Phase_Ins*(Mf<f1) + (Mf>=f1)*Phase_Int*(Mf<f2) + (Mf>=f2)*Phase_MR
	
	return Phase_IMR

@jit
def PhsD(eta,Coeffs,varphi,fD_QNM,rho_lm,tau_lm,Mf):
	f_RD, f_damp = fD_QNM

	f1 = 0.018
	f2 = 0.5*f_RD

	Sigma 	= Coeffs[7:11]
	Beta 	= Coeffs[11:14]
	Alpha 	= Coeffs[14:19]

	Connections = Phase_Connections(f1,f2,eta,varphi,Alpha,Beta,Sigma,f_RD,f_damp,rho_lm,tau_lm)
	Phase = Phase_IMR_D(eta,Coeffs,Connections,varphi,f1,f2,f_RD,f_damp,rho_lm,tau_lm,Mf)
	
	return Phase

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
def Amplitude_MR(x,eta,freq_QNM,gamma):
	f_RD, f_damp = freq_QNM
	gamma1, gamma2, gamma3 = gamma
	df = (x-f_RD)
	df2 = df*df
	arg = -gamma2*df/(gamma3*f_damp)
	K = funcA0(x,eta)
	A_MR = K*gamma1*gamma3*f_damp/(df2+(gamma3*f_damp)**2) * jnp.exp(arg)
	return A_MR

@jit
def diff_MR(x,eta,freq_QNM,gamma):
	f_RD, f_damp = freq_QNM
	gamma1, gamma2, gamma3 = gamma
	df = x-f_RD
	df2 = df*df
	D = df2+(gamma3*f_damp)**2
	Amr = Amplitude_MR(x,eta,freq_QNM,gamma)
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
def fpoly(ff, a):
	a = jnp.array(a,dtype=jnp.float64)
	n = jnp.arange(5)
	return a[:5] * ff ** n

@jit
def Amplitude_IMRPhenomD(*args_amp):
	dL, M, eta, chi_a, chi_s, deltaM, Coeffs, freq_QNM, Mf = args_amp

	eta2 = eta*eta
	eta3 = eta2*eta
	eta4 = eta3*eta

	f_RD, f_damp = freq_QNM

	#================================================
	# Coefficients
	#================================================

	# From Eqs. (B14) - (B20) of [PhysRevD.93.044007]

	Ai = [0]*7
	Ai[0] = 1.
	Ai[1] = 0.
	Ai[2] = -(323./224) + (451./168)*eta
	Ai[3] = 27*deltaM*chi_a/8 + (27./8 - 11*eta/6)*chi_s
	Ai[4] = -27312085./8128512 - 1975055*eta/338688 + 105271*eta2/24192 +\
	 (-81./32 + 8*eta)*chi_a**2 - 81.*deltaM*chi_a*chi_s/16 + (-81./32 + 17*eta/8)*chi_s**2
	Ai[5] = -85*PI/64 + 85*PI*eta/16 + deltaM*(285197./16128 - 1579*eta/4032)*chi_a \
	+ (285197./16128 - 15317*eta/672 - 2227*eta2/1008)*chi_s
	Ai[6] = -177520268561./8583708672 + (545384828789./5007163392 - 205*pi2/48)*eta \
	- 3248849057.*eta2/178827264 + 34473079*eta3/6386688 \
	+ (1614569./64512 - 1873643*eta/16128 + 2167*eta2/42)*chi_a**2 \
	+ (31*PI/12 - 7*PI*eta/3)*chi_s + (1614569./64512 - 61391*eta/1344 + 57451*eta2/4032)*chi_s**2 \
	+ deltaM*chi_a*(31*PI/12 + (1614569./32256 - 165961*eta/2688)*chi_s)
	
	gamma1, gamma2, gamma3 = Coeffs[4:7]
	rho1, rho2, rho3 = Coeffs[:3]
	u2 = Coeffs[3]

	rho = [rho1,rho2,rho3]
	gamma = [gamma1, gamma2, gamma3]
	
	f1 = 0.014
	f3 = jnp.abs( f_RD + f_damp*gamma3*(jnp.sqrt(1.-gamma2**2)-1.)/gamma2 )
	f2 = (f1+f3)/2

	u1 = Amplitude_Ins(f1,eta,Ai,rho)
	u3 = Amplitude_MR(f3,eta,freq_QNM,gamma)
	d1 = diff_Ins(f1,eta,Ai,rho)
	d3 = diff_MR(f3,eta,freq_QNM,gamma)

	K1 = funcA0(f1,eta)
	K2 = funcA0(f2,eta)
	K3 = funcA0(f3,eta)
	u2 *= K2

	r1 = fpoly(f1,jnp.ones(5))*K1
	r2 = fpoly(f2,jnp.ones(5))*K2
	r3 = fpoly(f3,jnp.ones(5))*K3
	r4 = fpoly(f1,[-7.,-1.,5.,11.,17.])*K1 / (6.*f1)
	r5 = fpoly(f3,[-7.,-1.,5.,11.,17.])*K3 / (6.*f3)

	Matrix = jnp.array([r1,r2,r3,r4,r5])
	delta = jnp.linalg.solve(Matrix,jnp.array([u1,u2,u3,d1,d3]))
	
	A_Ins = Amplitude_Ins(Mf,eta,Ai,rho)
	A_MR  = Amplitude_MR(Mf,eta,freq_QNM,gamma)
	A_Int = Amplitude_Int(Mf,eta,delta)

	A_IMR = A_Ins*(Mf<f1) + (Mf>=f1)*A_Int*(Mf<f3) + (Mf>=f3)*A_MR

	Amp = A0 * A_IMR * M*M*(GMc3**(7./6))/dL

	return Amp

#=====================================================#=====================================================
#----------------------------------------------Amplitude HM------------------------------------------------
#=====================================================#=====================================================

@jit
def func_Hlm(eta,chi_a,chi_s,deltaM,Mf,l,m):
	
	vm = jnp.array([(2.*PI*Mf/m)**(n/3) for n in range(1,5)],dtype=jnp.float64)
	v, v2, v3, v4 = vm
	args = [vm,deltaM,chi_a,chi_s,eta]

	# Define the branch functions
	branches = [
		lambda _: (Mf*0+1.).astype(jnp.complex128), # H22
		lambda _: (jnp.sqrt(2.)/3.*(
						(v*deltaM - 1.5*v2*(chi_a + deltaM*chi_s)+ v3*deltaM*(335./672 + 117./56*eta))\
						+ v4*(chi_a*(3427./1344 - 2101./336*eta) \
							+ deltaM*chi_s*(3427./1344 - 965./336*eta) \
							+ deltaM*(-0.5j - PI - 2.j * 0.69314718056) ) \
					) ).astype(jnp.complex128), # H21
		lambda _: (0.75 * jnp.sqrt(5./7.) * v * deltaM).astype(jnp.complex128), 			 # H33
		lambda _: (1./3. * jnp.sqrt(5./7.) * v2 * (1.-3.*eta)).astype(jnp.complex128), 		 # H32
		lambda _: (4./9. * jnp.sqrt(10./7.) * v2 * (1.-3.*eta)).astype(jnp.complex128), 	 # H44
		lambda _: (0.75 * jnp.sqrt(3./35.) * v3 * deltaM*(1.-2.*eta)).astype(jnp.complex128) # H43
	]
	lm_to_index = jnp.array([2+2,2+1,3+3,3+2,4+4,4+3])
	index = jnp.where(lm_to_index==(l+m),size=1)[0][0]
	branches.append(lambda args: 0.j*Mf)  # Default branch
	Hlm = jax.lax.switch(index, branches, args)
	
	return jnp.real( jnp.sqrt(Hlm*jnp.conj(Hlm)) )

@jit
def pow_kappaN(kp):
	return [kp**n for n in range(1,7)]

@jit
def funcZ_lm(m1,m2,chi_final):

	arg_kappa = jnp.log(2.0-chi_final) / jnp.log(3.)

	kappa = jnp.zeros([5,5])

	ll = jnp.array(idxs_HM)[:,0]
	mm = jnp.array(idxs_HM)[:,1]
	@jit
	def compute_Klm(ll,mm):
		return arg_kappa ** (1./(2.+ll-jnp.abs(mm)))
	
	vectorized_compute_Klm = vmap(compute_Klm,in_axes=(0,0))
	@jit
	def update_Klm(kappa):
		kappa_val = vectorized_compute_Klm(ll,mm)
		return kappa.at[ll,mm].set(kappa_val)
	kappa = update_Klm(kappa)

	Z = jnp.zeros([5,5],dtype=complex)
	# non-null idxs_HM = [[2,2],[2,1],[3,3],[3,2],[4,4],[4,3]]

	l,m = 2,2 ; kp1, kp2, kp3, kp4, kp5, kp6 = [kappa[l][m]**n for n in range(1,7)]
	Zlm = 1. + \
			kp1 * 1.557847   * jnp.exp(2.903124j) + \
			kp2 * 1.95097051 * jnp.exp(5.920970j) + \
			kp3 * 2.09971716 * jnp.exp(2.760585j) + \
			kp4 * 1.41094660 * jnp.exp(5.914340j) + \
			kp5 * 0.41063923 * jnp.exp(2.795235j)
	Z = Z.at[l,m].set(Zlm)

	l,m = 3,2 ; kp1, kp2, kp3, kp4, kp5, kp6 = [kappa[l][m]**n for n in range(1,7)]
	Zlm = 1.022464*jnp.exp(0.004870j) + \
			kp1 * 0.24731213  * jnp.exp(0.665292j) + \
			kp2 * 1.70468239  * jnp.exp(3.138283j) + \
			kp3 * 0.94604882  * jnp.exp(0.163247j) + \
			kp4 * 1.53189884  * jnp.exp(5.703573j) + \
			kp5 * 2.28052668  * jnp.exp(2.685231j) + \
			kp6 * 0.92150314  * jnp.exp(5.841704j)
	Z = Z.at[l,m].set(Zlm)

	l,m = 4,4 ; kp1, kp2, kp3, kp4, kp5, kp6 = [kappa[l][m]**n for n in range(1,7)]
	Zlm = 2. + \
			kp1 * 2.658908   * jnp.exp(3.002787j) + \
			kp2 * 2.97825567 * jnp.exp(6.050955j) + \
			kp3 * 3.21842350 * jnp.exp(2.877514j) + \
			kp4 * 2.12764967 * jnp.exp(5.989669j) + \
			kp5 * 0.60338186 * jnp.exp(2.830031j)
	Z = Z.at[l,m].set(Zlm)

	l,m = 2,1 ; kp1, kp2, kp3, kp4, kp5, kp6 = [kappa[l][m]**n for n in range(1,7)]
	Zlm = 0.589113 * jnp.exp(0.043525j) + \
			kp1 * 0.18896353  * jnp.exp(2.289868j) + \
			kp2 * 1.15012965  * jnp.exp(5.810057j) + \
			kp3 * 6.04585476  * jnp.exp(2.741967j) + \
			kp4 * 11.12627777 * jnp.exp(5.844130j) + \
			kp5 * 9.34711461  * jnp.exp(2.669372j) + \
			kp6 * 3.03838318  * jnp.exp(5.791518j)
	Z = Z.at[l,m].set(Zlm)

	l,m = 3,3 ; kp1, kp2, kp3, kp4, kp5, kp6 = [kappa[l][m]**n for n in range(1,7)]
	Zlm = 1.5 + \
			kp1 * 2.095657   * jnp.exp(2.964973j) + \
			kp2 * 2.46964352 * jnp.exp(5.996734j) + \
			kp3 * 2.66552551 * jnp.exp(2.817591j) + \
			kp4 * 1.75836443 * jnp.exp(5.932693j) + \
			kp5 * 0.49905688 * jnp.exp(2.781658j)
	Z = Z.at[l,m].set(Zlm)

	l,m = 4,3 ; kp1, kp2, kp3, kp4, kp5, kp6 = [kappa[l][m]**n for n in range(1,7)]
	Zlm = 1.5 + \
			kp1 * 0.205046   * jnp.exp(0.595328j) + \
			kp2 * 3.10333396 * jnp.exp(3.016200j) + \
			kp3 * 4.23612166 * jnp.exp(6.038842j) + \
			kp4 * 3.02890198 * jnp.exp(2.826239j) + \
			kp5 * 0.90843949 * jnp.exp(5.915164j)
	Z = Z.at[l,m].set(Zlm)

	return Z

@jit
def freq_RD_Damp_lm(m1,m2,chi_final,M_final):
	Zlm     = funcZ_lm(m1,m2,chi_final)
	inv2PiM = 1./(2*PI*M_final)
	
	return inv2PiM * Zlm.real, inv2PiM * Zlm.imag

@jit
def func_fAmp(f_ring_lm,Mf):
	ones = jnp.ones([5,5])
	zeros = jnp.zeros([5,5])
	
	f_ring_22 = f_ring_lm[2][2]
	rho_lm    = f_ring_22 / f_ring_lm 

	f0 = 0.014

	f1 = f0*ones
	fi = f1/rho_lm
	fr = f_ring_lm

	try:
		Empty = jnp.zeros([5,5,len(Mf)], dtype=jnp.float64)
	except:
		Empty = jnp.zeros([5,5], dtype=jnp.float64)

	# idxs_HM = [[2,2],[2,1],[3,3],[3,2],[4,4],[4,3]]
	ll = jnp.array(idxs_HM)[:,0]
	mm = jnp.array(idxs_HM)[:,1]
	
	@jit
	def compute_famp(ll,mm):
		fr_lm = fr[ll][mm]
		fi_lm = f0*fr_lm/f_ring_22
		
		Ci = 2.*Mf/mm
		Cm = (f_ring_22-2.*fi_lm/mm)/(fr_lm-fi_lm) * (Mf-fi_lm) + 2.*fi_lm/mm
		Cr = Mf - (fr_lm - f_ring_22)

		return Ci*(Mf<fi_lm) + (fi_lm<=Mf)*Cm*(Mf<fr_lm) + (fr_lm<=Mf)*Cr
	
	vectorized_compute_famp = vmap( compute_famp, in_axes=(0,0) )

	@jit
	def update_famp(f_amp):
		val = vectorized_compute_famp(ll,mm)
		return f_amp.at[ll,mm].set(val)
	f_amp = update_famp(Empty)

	return f_amp

@jit
def AmpHM_lm(Mc,eta,dL,chi_a,chi_s,Coeffs,f_QNM,Mf):
	M = Mc/eta**(3./5)
	deltaM = jnp.sqrt(1. - 4*eta)
	
	_, _, f_ring_lm, _ = f_QNM
	
	f_amp   = func_fAmp(f_ring_lm,Mf)
	try:
		empty_complex = jnp.zeros([5,5,len(Mf)],dtype=complex) # Hlm_term
	except:
		empty_complex = jnp.zeros([5,5],dtype=complex) # Hlm_term
	f_amp_22 = f_amp[2][2]
	
	ll = jnp.array(idxs_HM)[:,0]
	mm = jnp.array(idxs_HM)[:,1]

	@jit
	def compute_Alm(l,m):
		coeff = 1.

		H22_fa  = func_Hlm(eta,0,0,deltaM,f_amp[l][m],2,2) # := 1
		Hlm_fa  = func_Hlm(eta,chi_a,chi_s,deltaM,f_amp[l][m],l,m)
		Hlm_Mf  = func_Hlm(eta,chi_a,chi_s,deltaM,Mf,l,m)
		Hlm_2Mf = func_Hlm(eta,chi_a,chi_s,deltaM,2*Mf/m,l,m)

		return jnp.where(Hlm_Mf == 0, Hlm_Mf, coeff * Hlm_fa * Hlm_Mf / Hlm_2Mf / H22_fa) # Avoid problems when m1=m2 and (l,m) in [(3,3),(4,3)]

	vectorized_compute_Alm = vmap( compute_Alm, in_axes=(0,0) )
	@jit
	def update_Alm(coeff_amp):
		Alm = vectorized_compute_Alm(ll,mm)
		return coeff_amp.at[ll,mm].set(Alm)
	coeff_amp = empty_complex.copy()
	coeff_amp = update_Alm(coeff_amp)
	#******************************************************
	# Aternative way to compute Amplitude(D)
	# Here we compute only the choosen (l,m) components
	# instead all (5x5) components
	#******************************************************
	f_RD, f_damp, f_ring_lm, f_damp_lm = f_QNM
	freq_QNM = [f_RD, f_damp]
	@jit
	def compute_AmpD(l,m):
		return Amplitude_IMRPhenomD(dL,M,eta,chi_a,chi_s,deltaM,Coeffs,freq_QNM,f_amp[l][m])

	vectorized_compute_AmpD = vmap( compute_AmpD, in_axes=(0,0) )
	@jit
	def update_AmpD(tmp):
		AmpD_vals = vectorized_compute_AmpD(ll,mm)
		return tmp.at[ll,mm].set( AmpD_vals )
	AmpD = empty_complex.copy()
	AmpD = update_AmpD(AmpD)
	
	return coeff_amp * AmpD

#=====================================================#=====================================================
#----------------------------------------------Phase HM------------------------------------------------
#=====================================================#=====================================================

@jit
def fAB(a,b,Mf):
	return a*jnp.array(Mf)+b

@jit
def PhaseD_lm(eta,varphi,Coeffs,f_qnm,rho_lm,tau_lm,ll,mm,Mf):

	f_RD, f_damp, f_ring_lm, f_ring_22 = f_qnm
	fD_QNM = [f_RD, f_damp]

	f0phs = 0.018

	f1 = f0phs
	fi_lm = f1/rho_lm
	fr_lm = f_ring_lm
	
	args_lm = [eta,Coeffs,varphi,fD_QNM,rho_lm,tau_lm]

	#**********************************************************************************
	ai, bi = 2./mm, 0
	am = (f_ring_22-2*fi_lm/mm)/(fr_lm-fi_lm)
	bm = fi_lm*(2./mm - am)
	ar, br = rho_lm, 0

	Mfi = fAB(ai,bi,Mf); Mfii = fAB(ai,bi,fi_lm)
	Mfm = fAB(am,bm,Mf); Mfmi = fAB(am,bm,fi_lm); Mfmr = fAB(am,bm,fr_lm)
	Mfr = fAB(ar,br,Mf); Mfrr = fAB(ar,br,fr_lm)

	C1_lm = PhsD(*args_lm,Mfii)/ai - PhsD(*args_lm,Mfmi) /am
	C2_lm = PhsD(*args_lm,Mfmr)/am - PhsD(*args_lm,Mfrr)/ar + C1_lm

	# PhsD(m1,m2,Coeffs,varphi,fD_QNM,rho_lm,tau_lm,Mf)
	Phase_ins = PhsD(*args_lm,Mfi)/ai
	Phase_int = PhsD(*args_lm,Mfm)/am + C1_lm
	Phase_mr  = PhsD(*args_lm,Mfr)/ar + C2_lm

	cShift = jnp.array([0,PI/2,0,-PI/2,PI],dtype=jnp.float64)
	phase_lm = cShift[mm] + Phase_ins*(Mf<fi_lm) + (fi_lm<=Mf)*Phase_int*(Mf<fr_lm) + (fr_lm<=Mf)*Phase_mr
	#**********************************************************************************
	#phase_lm = PhsD(*args_lm, Mf)
	
	return phase_lm

# PhsD(mass1,mass2,sz1,sz2,t0,phi0,fRef,rho_lm,tau_lm,freq)
@jit
def PhaseHM_lm(eta,deltaM,chi_a,chi_s,Mf_Ref,Coeffs,f_QNM,Rho,Tau,Mf):
 
	gamma2, gamma3 = Coeffs[5:7]
	Sigma = Coeffs[7:11]
	Beta  = Coeffs[12:14]
	Alpha = Coeffs[14:19]

	alpha5 = Alpha[4]
	gamma22 = gamma2*gamma2

	f_RD, f_damp, f_ring_lm, f_damp_lm = f_QNM
	
	f_peak = jnp.abs( f_RD + (f_damp*gamma3/gamma2 )*(jnp.sqrt(1.-gamma22)-1.))
	#-------------------------------------------------------------

	T0 = Diff_PhiD_MR(f_peak,eta,Alpha,f_RD,f_damp,1,1) # rho_lm = tau_lm = 1
	phi1 = -(Mf - Mf_Ref)*T0

	varphi = VarPhi_InsPhase(eta,deltaM,chi_a,chi_s)

	#*********************************************************************************
	f_qnm_22 = f_RD, f_damp, f_ring_lm[2][2], f_ring_lm[2][2]
	pH0 = PhaseD_lm(eta,varphi,Coeffs,f_qnm_22,1.,1.,2,2,Mf_Ref)/2
	@jit
	def Compute_PhaseHM(l,m):
		rho_lm = Rho[l][m]
		tau_lm = Tau[l][m]
		f_qnm = f_RD, f_damp, f_ring_lm[l][m], f_ring_lm[2][2]
		# PhaseD_lm(m1,m2,varphi,Coeffs,f_qnm,rho_lm,tau_lm,ll,mm,Mf)

		#pH0 = PhaseD_lm(m1,m2,varphi,Coeffs,f_qnm,1.,1.,2,2,Mf_Ref)/2
		phase_val = phi1 - (m * pH0) + PhaseD_lm(eta,varphi,Coeffs,f_qnm,rho_lm,tau_lm,l,m,Mf)
		return phase_val

	vectorized_compute_PhaseHM = vmap( Compute_PhaseHM, in_axes=(0,0) )
	
	ll = jnp.array(idxs_HM)[:,0]
	mm = jnp.array(idxs_HM)[:,1]

	PhaseHM_vals = vectorized_compute_PhaseHM( ll, mm )

	try:
		PhaseHM = jnp.zeros([5,5,len(Mf)],dtype=jnp.float64)
	except:
		PhaseHM = jnp.zeros([5,5],dtype=jnp.float64)

	PhaseHM = PhaseHM.at[ll,mm].set(PhaseHM_vals) 
	#*********************************************************************************
	return PhaseHM 

@jit
def SphericalHarm(theta,phi,l,m):
	branches = [
		lambda _:    jnp.sqrt(5./(64*PI))*(1+jnp.cos(theta))**2 ,
		lambda _:    jnp.sqrt(5./(16*PI))*jnp.sin(theta)*(1+jnp.cos(theta)) ,
		lambda _:    jnp.sqrt(5./(16*PI))*jnp.sin(theta)*(1-jnp.cos(theta)) ,
		lambda _:    jnp.sqrt(5./(64*PI))*(1-jnp.cos(theta))**2 ,
		lambda _:   -jnp.sqrt(21./(2*PI))*jnp.cos(theta/2)**5*jnp.sin(theta/2) ,
		lambda _:    jnp.sqrt(7./PI)*jnp.cos(theta/2)**4*(3*jnp.cos(theta)-2)/2 ,
		lambda _:    jnp.sqrt(7./PI)*jnp.sin(theta/2)**4*(3*jnp.cos(theta)+2)/2 ,
		lambda _:    jnp.sqrt(21./(2*PI))*jnp.sin(theta/2)**5*jnp.cos(theta/2) ,
		lambda _:  3*jnp.sqrt(7./PI)*jnp.cos(theta/2)**6*jnp.sin(theta/2)**2 ,
		lambda _: -3*jnp.sqrt(7/(2*PI))*jnp.cos(theta/2)**5*(2*jnp.cos(theta)-1)*jnp.sin(theta/2) ,
		lambda _:  3*jnp.sqrt(7/(2*PI))*jnp.sin(theta/2)**5*(2*jnp.cos(theta)+1)*jnp.cos(theta/2) ,
		lambda _:  3*jnp.sqrt(7./PI)*jnp.sin(theta/2)**6*jnp.cos(theta/2)**2
	]

	Y_lm_idxs = [[2,2],[2,1],[2,-1],[2,-2],[3,3],[3,2],[3,-2],[3,-3],[4,4],[4,3],[4,-3],[4,-4]]

	idxs = jnp.array( [ 10*x+y for x,y in Y_lm_idxs ] )
	index = jnp.where(idxs==(10*l+m),size=1)[0][0]
	branches.append(lambda _: 0.)  # Default branch
	Ylm = jax.lax.switch(index, branches, 0)

	Ylm *= jnp.exp(1.j*phi*m)
	return Ylm

@jit
def hphx_IMRPhenomHM(dL,iota,phi0,t0,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,freq): # 12 GwPrms + freq
	chi1, chi2 = sz1, sz2
	m1 = 0.5*(Mc/eta**(3./5)) * (1 + jnp.sqrt(1-4*eta))
	m2 = 0.5*(Mc/eta**(3./5)) * (1 - jnp.sqrt(1-4*eta))
	M = m1+m2
	deltaM = (m1-m2)/M

	GM = M*GMc3
	Mf = freq*GM
	fRef = 1.0
	Mf_Ref = GM*fRef

	ll = jnp.array(idxs_HM)[:,0]
	mm = jnp.array(idxs_HM)[:,1]

	#================================================
	# Ringdown and Damping Frequencies
	#================================================
	
	S = (m1*m1*chi1 + m2*m2*chi2)/M**2
	S_hat = S_hat = (m1*m1*chi1 + m2*m2*chi2)/(m1*m1+m2*m2)

	#af, E_rad = Final_Spin(eta,S,S_hat)
	eta2 = eta*eta
	eta3 = eta2*eta
	eta4 = eta3*eta

	#================================================
	# Final Spin: Eq. 3.6 of arXiv:1508.07250
	#================================================
	S2 = S*S
	S3 = S2*S
	S4 = S3*S

	chi_final = eta*(3.4641016151377544 - 4.399247300629289*eta + \
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
	
	M_final = 1.- E_rad

	f_RD   = func_interp_RD(chi_final) / M_final  
	f_damp = func_interp_damp(chi_final) / M_final 

	#f_RD   = JAX_CubicSpline(af, aa, fRg) / (1.-E_rad) # Recovery lalsimulation
	#f_damp = JAX_CubicSpline(af, aa, fDp) / (1.-E_rad) # Recovery lalsimulation

	#================================================

	#f_ring_lm, f_damp_lm = freq_RD_Damp_lm(m1,m2,chi1,chi2) 
	f_ring_lm, f_damp_lm = freq_RD_Damp_lm(m1,m2,chi_final,M_final)
	f_QNM = [f_RD, f_damp, f_ring_lm, f_damp_lm]

	f_ring_22 = f_ring_lm[2][2]
	f_damp_22 = f_damp_lm[2][2]

	#Rho  = f_ring_22 / f_ring_lm
	#Tau  = f_damp_lm / f_damp_22

	Empty = jnp.zeros([5,5])
	Rho = Empty.copy()
	Tau = Empty.copy()

	def Compute_Rho(l,m):
		return f_ring_22 / f_ring_lm[l][m]
	def Compute_Tau(l,m):
		return f_damp_lm[l][m] / f_damp_22

	vectorized_Rho = vmap(Compute_Rho)
	vectorized_Tau = vmap(Compute_Tau)
	vals_Rho = vectorized_Rho(ll,mm)
	vals_Tau = vectorized_Tau(ll,mm)

	Rho = Rho.at[ll,mm].set(vals_Rho)
	Tau = Tau.at[ll,mm].set(vals_Tau)

	#-------------------------------------------------------------
	# Computing Alpha,Beta,Gamma,Sigma Coefficients:
	#-------------------------------------------------------------
	xi_eff = (chi1*m1 + chi2*m2)/M
	xi_PN = xi_eff - (38./113.)*eta*(chi1+chi2)

	dXi = xi_PN - 1.
	Coeffs = vmap( lambda i: Lambda(i, eta, dXi) )(jnp.arange(19))
	# rho's 	= Coeffs[:3]
	# v2 		= Coeffs[3]
	# gamma's 	= Coeffs[4:7]
	# sigma's 	= Coeffs[7:11]
	# beta's 	= Coeffs[11:14]
	# alpha's 	= Coeffs[14:]

	#-------------------------------------------------------------

	empty_lm = f_ring_lm*0
	ones_lm = empty_lm + 1
	
	chi_a = 0.5*(chi1 - chi2)
	chi_s = 0.5*(chi1 + chi2)

	AmpHM = AmpHM_lm(Mc,eta,dL,chi_a,chi_s,Coeffs,f_QNM,Mf)
	PhsHM = PhaseHM_lm(eta,deltaM,chi_a,chi_s,Mf_Ref,Coeffs,f_QNM,Rho,Tau,Mf) + (2*PI*freq*t0 - phi0) # Check carrefully this last term (2*PI*freq*t0 - phi0)

	# idxs_HM = [[2,2],[2,1],[3,3],[3,2],[4,4],[4,3]]
	k = 2*jnp.sqrt(5./(64*PI))
	@jit
	def compute_hlm(l,m):
		Ylm1 = SphericalHarm(iota,0,l,m)
		Ylm2 = jnp.conj(SphericalHarm(iota,0,l,-m))

		gp_lm =  0.5  * ( Ylm1 + (-1)**l * Ylm2) / k
		gx_lm =  -0.5j * ( Ylm1 - (-1)**l * Ylm2) / k

		Phi_lm = PhsHM[l][m]
		Amp_lm = AmpHM[l][m]
		
		hlm = Amp_lm * jnp.exp( -1.j * Phi_lm )
		
		Hp_lm = gp_lm * hlm
		Hx_lm = gx_lm * hlm
		return Hp_lm, Hx_lm
	vectorized_compute_hlm = vmap(compute_hlm, in_axes=(0,0))
	@jit
	def update_hp_hx(ll,mm):
		try:
			hp = jnp.zeros(len(freq),dtype=complex)
			hx = jnp.zeros(len(freq),dtype=complex)
		except:
			hp, hx = 0.j, 0.j
		Hp_lm_values, Hx_lm_values = vectorized_compute_hlm(ll,mm)
		hp = jnp.sum(Hp_lm_values, axis=0)
		hx = jnp.sum(Hx_lm_values, axis=0)
		return hp, hx
	hp, hx = update_hp_hx(ll,mm)

	return hp, hx
