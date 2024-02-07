import numpy as np
import GWDALI.lib.Waveforms as wf
import GWDALI.lib.Angles_lib as geo
import GWDALI.lib.Dictionaries as gwdict
from scipy.interpolate import interp1d

rad = np.pi/180
deg = 1./rad

c = 299792458 # m/s
R_earth = 6371.e3 # meters

PSD, labels_tex = gwdict.Load_Dictionaries()


def Pattern_Func(alpha,beta,psi,Omega):
	u = np.cos(beta) ; Coeff = np.sin(Omega)
	fp = 0.5*(1.+u**2)*np.cos(2.*alpha)*Coeff
	fx = -u*np.sin(2.*alpha)*Coeff
	
	Fp =  fp*np.cos(2.*psi) + fx*np.sin(2.*psi)
	Fx = -fp*np.sin(2.*psi) + fx*np.cos(2.*psi)

	return Fp, Fx

#-------------------------------------------------
def Integral(x,y):
	return np.sum( np.array( [0.5*(y[i]+y[i-1])*(x[i]-x[i-1]) for i in range(1,len(x)) ] ) )

def ScalarProduct(freq,Sn,A,B):
	return 4*np.real( Integral(freq, A*np.conj(B)/Sn) )

def GW_Polarizations(params, freq, approx):
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

	hp, hx, freq0 = wf.Waveforms(m1,m2,iota,DL,s1,s2,freq, approx=approx)

	return hp, hx, freq0

def Signal(params, det, approx):
	alpha    = params['RA']*rad       # rad
	beta     = (90-params['Dec'])*rad # rad
	iota     = params['iota']         # rad
	psi      = params['psi']          # rad
	t_coal   = params['t_coal']       # sec
	phi_coal = params['phi_coal']     # rad
	
	name  = det['name']
	freq  = det['freq'].copy()
	lon   = det['lon'] # deg
	lat   = det['lat'] # deg
	rot   = det['rot'] # deg

	alpha_obs, beta_obs, psi_obs = geo.ObsAngles(alpha,beta,iota,psi,lon,lat,rot)

	t_delay = -np.cos(beta_obs)*R_earth/c # time delay between detector and center of Earth

	Fp, Fx = Pattern_Func(alpha_obs,beta_obs,psi_obs,det['shape']*rad)
	hp, hx, freq0 = GW_Polarizations(params, freq, approx)

	Phase = 2*np.pi*freq0*(t_coal+t_delay) - phi_coal
	H = (Fp*hp + Fx*hx)*np.exp(1.j*Phase)
	
	gw_signal = interp1d(freq0,H,bounds_error=False,fill_value='extrapolate')
	h = gw_signal(freq)

	return h	

#-------------------------------------------------

# eps (standard) = 1.e-6
def split_prms(params,x,eps, diff_order):
	p = params[x]
	dx = np.max([eps,eps*p])
	P0 = params.copy() ; P0[x] = p - dx	
	P1 = params.copy() ; P1[x] = p - dx/2
	P2 = params.copy() ; P2[x] = p + dx/2
	P3 = params.copy() ; P3[x] = p + dx
	if(diff_order == 2): return [P1,P2], dx
	elif(diff_order==4): return [P0,P1,P2,P3], dx
	else:
		print("\n\t Invalid diff_order! Allowed values: [2,4] ")
		quit()

def Diff1(x, params, det, approx, eps, diff_order):
	Ps, dx = split_prms(params,x,eps, diff_order)
	Y = [Signal(P, det, approx) for P in Ps]
	if(diff_order == 2): return (Y[1]-Y[0])/dx
	elif(diff_order==4): return 4*(Y[2]-Y[1])/(3*dx) - (Y[3]-Y[0])/(6*dx)

def Diff2(xi, xj, params, det, approx, eps, diff_order):
	Ps, dx = split_prms(params,xi,eps, diff_order)
	Y = [Diff1(xj, P, det, approx, eps, diff_order) for P in Ps]
	if(diff_order == 2): return (Y[1]-Y[0])/dx
	elif(diff_order==4): return 4*(Y[2]-Y[1])/(3*dx) - (Y[3]-Y[0])/(6*dx)

def Diff3(xi, xj, xk, params, det, approx, eps, diff_order):
	Ps, dx = split_prms(params,xi,eps, diff_order)
	Y = [Diff2(xj,xk, P, det, approx, eps, diff_order) for P in Ps]
	if(diff_order == 2): return (Y[1]-Y[0])/dx
	elif(diff_order==4): return 4*(Y[2]-Y[1])/(3*dx) - (Y[3]-Y[0])/(6*dx)

#-------------------------------------------------#-------------------------------------------------
#-------------------------------------------------#-------------------------------------------------

def Fisher_ij(xi,xj, params, det, approx, eps, diff_order): # [1,1]
	Dxi = Diff1(xi, params, det, approx, eps, diff_order)
	Dxj = Diff1(xj, params, det, approx, eps, diff_order)
	return ScalarProduct(det['freq'], det['Sn'], Dxi, Dxj )

#-------------------------------------------------#-------------------------------------------------
# (arXiv:2203.02670)

def Doublet3(xi,xj,xk, params, det, approx, eps, diff_order): # [1,2]
	D_i  = Diff1(xi, params, det, approx, eps, diff_order)
	D_jk = Diff2(xj,xk, params, det, approx, eps, diff_order)
	return ScalarProduct(det['freq'], det['Sn'], D_i, D_jk)

def Doublet4(xi,xj,xk,xl, params, det, approx, eps, diff_order): # [2,2]
	D_ij = Diff2(xi,xj, params, det, approx, eps, diff_order)
	D_kl = Diff2(xk,xl, params, det, approx, eps, diff_order)
	return ScalarProduct(det['freq'], det['Sn'], D_ij, D_kl)

#-------------------------------------------------#-------------------------------------------------

def Triplet4(xi,xj,xk,xl, params, det, approx, eps, diff_order): # [1,3]
	D_i   = Diff1(xi, params, det, approx, eps, diff_order)
	D_jkl = Diff3(xj,xk,xl, params, det, approx, eps, diff_order)
	return ScalarProduct(det['freq'], det['Sn'], D_i, D_jkl)

def Triplet5(xi,xj,xk,xl,xm, params, det, approx, eps, diff_order): # [2,3]
	D_ij  = Diff2(xi,xj, params, det, approx, eps, diff_order)
	D_klm = Diff3(xk,xl,xm, params, det, approx, eps, diff_order)
	return ScalarProduct(det['freq'], det['Sn'], D_ij, D_klm)

def Triplet6(xi,xj,xk,xl,xm,xn, params, det, approx, eps, diff_order): # [3,3]
	D_ijk = Diff3(xi,xj,xk, params, det, approx, eps, diff_order)
	D_lmn = Diff3(xl,xm,xn, params, det, approx, eps, diff_order)
	return ScalarProduct(det['freq'], det['Sn'], D_ijk, D_lmn)

#-------------------------------------------------#-------------------------------------------------
