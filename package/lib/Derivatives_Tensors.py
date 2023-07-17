import numpy as np
import GWDALI.lib.Waveforms as wf
import GWDALI.lib.Angles_lib as geo
import GWDALI.lib.Dictionaries as gwdict

rad = np.pi/180
deg = 1./rad

c = 3.e8
R_earth = 6371.e3

PSD, labels_tex = gwdict.Load_Dictionaries()

def Pattern_Func(alpha,beta,psi,shape): # signs fixed
	Coeff = np.sin(shape) #np.sqrt(3.)/2
	fp = Coeff*( 0.5*(1.+np.cos(beta)**2)*np.cos(2.*alpha) )
	fx = -Coeff*( np.cos(beta)*np.sin(2.*alpha) )
	
	Fp =  fp*np.cos(2.*psi) + fx*np.sin(2.*psi)
	Fx = -fp*np.sin(2.*psi) + fx*np.cos(2.*psi)
	#print('alpha,beta,psi:%.2f, %.2f, %.2f'% (alpha,beta,psi), end='\t')
	#print('\nFp, Fx:', Fp, Fx,'\tFp^2+Fx^2',Fp**2+Fx**2,'\t psi_det:',psi,'\n')
	#print('Fp^2+Fx^2:%.2f' % (Fp**2+Fx**2) , ' ; Fp, Fx: %.2f, %.2f' %(Fp, Fx))
	return Fp, Fx

#-------------------------------------------------
def Integral(x,y):
	return np.sum( np.array( [0.5*(y[i]+y[i-1])*(x[i]-x[i-1]) for i in range(1,len(x)) ] ) )

def ScalarProduct(freq,Sn,A,B):
	return 4*np.real( Integral(freq, A*np.conj(B)/Sn) )

def GW_Polarizations(params, freq, approximant):
	keys = list(params.keys())
	
	#-------------------------------------------------
	# Convert Masses
	#-------------------------------------------------
	if(all( p in  keys for p in ['m1','m2'])):
		m1 = params['m1'] ; m2 = params['m2']
	elif(all( p in keys for p in ['Mc','eta'] )):
		Mc = params['Mc'] ; eta = params['eta'] 
		Coeff = Mc/(2*eta**0.6)
		m1 = Coeff*( 1. + np.sqrt(1.-4*eta) ) # m1>m2
		m2 = Coeff*( 1. - np.sqrt(1.-4*eta) )
	elif(all( p in keys for p in ['Mc','q'] )):
		Mc = params['Mc'] ; q = params['q'] # q = m2/m1
		m2 = ( q*q*(1.+q) )**0.2 * Mc
		m1 = m2/q 
	else: print("\n\n# Invalid options to mass! \n\t >> Available: [m1,m2] or [Mc,eta] or [Mc,q]\n\n")
	#-------------------------------------------------

	DL       = params['DL']*1.e3      # Gpc --> Mpc
	iota     = params['iota']         # rad
	psi      = params['psi']          # rad

	sx1 = params['sx1']
	sy1 = params['sy1']
	sz1 = params['sz1']
	sx2 = params['sx2']
	sy2 = params['sy2']
	sz2 = params['sz2']

	s1 = [sx1, sy1, sz1]
	s2 = [sx2, sy2, sz2]

	hp, hx, _ = wf.Waveforms(m1,m2,iota,DL,s1,s2,freq, approx=approximant)

	return hp, hx

def Signal(params, detector, approximant):
	alpha    = params['RA']*rad       # rad
	beta     = (90-params['Dec'])*rad # rad
	iota     = params['iota']         # rad
	psi      = params['psi']          # rad
	t_coal   = params['t_coal']       # sec
	phi_coal = params['phi_coal']     # rad
	
	name  = detector['name']
	freq  = detector['freq'].copy()
	lon   = detector['lon'] # deg
	lat   = detector['lat'] # deg
	rot   = detector['rot']*rad # rad

	phi   = lon*rad
	theta = (90-lat)*rad

	alpha_obs, beta_obs = geo.AngTransf(alpha,beta,theta,phi,rot)
	psi_obs             = geo.poll_ang(alpha,beta,iota,psi,theta,phi,rot)
	#psi_obs = psi # Adding by hand!!!

	t_delay = -np.cos(beta_obs)*R_earth/c

	Fp, Fx = Pattern_Func(alpha_obs,beta_obs,psi_obs,detector['shape']*rad)
	hp, hx = GW_Polarizations(params, freq, approximant)

	exp2 = np.exp(-1.j*phi_coal)

	h = Fp*hp + Fx*hx
	
	return h * exp2	

#-------------------------------------------------

eps = 1.e-6
def split_prms(params,x):
	p = params[x]
	dx = np.max([1.e-10,eps*p])
	P1 = params.copy() ; P1[x] = p - dx/2
	P2 = params.copy() ; P2[x] = p + dx/2
	return P1, P2, dx

def Diff1(x,params , detector, approximant):
	P1, P2, dx = split_prms(params,x)
	y2 = Signal(P2, detector, approximant)
	y1 = Signal(P1, detector, approximant)
	return (y2-y1)/dx

def Diff2(xi,xj,params , detector, approximant):
	P1, P2, dx = split_prms(params,xi)
	y2 = Diff1(xj, P2, detector, approximant)
	y1 = Diff1(xj, P1, detector, approximant)
	return (y2-y1)/dx

def Diff3(xi,xj,xk,params , detector, approximant):
	P1, P2, dx = split_prms(params,xi)
	y2 = Diff2(xj,xk, P2, detector, approximant)
	y1 = Diff2(xj,xk, P1, detector, approximant)
	return (y2-y1)/dx

#-------------------------------------------------#-------------------------------------------------
#-------------------------------------------------#-------------------------------------------------

def Fisher_ij(xi,xj, params, detector, approximant): # [1,1]
	diff_xi = Diff1(xi,params , detector, approximant)
	diff_xj = Diff1(xj,params , detector, approximant)
	return ScalarProduct(detector['freq'], detector['Sn'], diff_xi, diff_xj )

#-------------------------------------------------#-------------------------------------------------

def func_doublet3(xi,xj,xk, params, detector, approximant): # [1,2]
	diff_xi    = Diff1(xi,params , detector, approximant)
	diff_xj_xk = Diff2(xj,xk,params , detector, approximant)
	return ScalarProduct(detector['freq'], detector['Sn'], diff_xi, diff_xj_xk)

def func_doublet4(xi,xj,xk,xl, params, detector, approximant): # [2,2]
	diff_xi_xj = Diff2(xi,xj,params , detector, approximant)
	diff_xk_xl = Diff2(xk,xl,params , detector, approximant)
	return ScalarProduct(detector['freq'], detector['Sn'], diff_xi_xj, diff_xk_xl)

#-------------------------------------------------#-------------------------------------------------

def func_triplet4(xi,xj,xk,xl, params, detector, approximant): # [1,3]
	diff_xi       = Diff1(xi,params , detector, approximant)
	diff_xj_xk_xl = Diff3(xj,xk,xl,params , detector, approximant)
	return ScalarProduct(detector['freq'], detector['Sn'], diff_xi, diff_xj_xk_xl)

def func_triplet5(xi,xj,xk,xl,xm, params, detector, approximant): # [2,3]
	diff_xi_xj    = Diff2(xi,xj,params , detector, approximant)
	diff_xk_xl_xm = Diff3(xk,xl,xm,params , detector, approximant)
	return ScalarProduct(detector['freq'], detector['Sn'], diff_xi_xj, diff_xk_xl_xm)

def func_triplet6(xi,xj,xk,xl,xm,xn, params, detector, approximant): # [3,3]
	diff_xi_xj_xk = Diff3(xi,xj,xk,params , detector, approximant)
	diff_xl_xm_xn = Diff3(xl,xm,xn,params , detector, approximant)
	return ScalarProduct(detector['freq'], detector['Sn'], diff_xi_xj_xk, diff_xl_xm_xn)

#-------------------------------------------------#-------------------------------------------------