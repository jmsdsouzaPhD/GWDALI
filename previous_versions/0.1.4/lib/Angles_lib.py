import numpy as np
import matplotlib.pyplot as plt

def scl_prod(V1,V2): # Scalar Product
	if(len(V1)!=len(V2)):
		print("Error!")
		quit()
	SP = 0 ; ndim = len(V1)
	for i in range(ndim): SP += V1[i]*V2[i]
	return SP

def vec_prod(V1,V2): # Vectorial Product
	if(len(V1)!=len(V2)):
		print("Error!")
		quit()
	vPx = V1[1]*V2[2] - V1[2]*V2[1]
	vPy = V1[2]*V2[0] - V1[0]*V2[2]
	vPz = V1[0]*V2[1] - V1[1]*V2[0]
	return vPx, vPy, vPz

def AngTransf(alpha,beta,theta,phi,ksi): # Sky Position from Geocentric to Detector Coordinates
	sA   = np.sin(beta)*np.sin(alpha-phi)
	cA   = np.sin(beta)*np.cos(theta)*np.cos(alpha-phi) - np.cos(beta)*np.sin(theta)
	cosB = np.sin(beta)*np.sin(theta)*np.cos(alpha-phi) + np.cos(beta)*np.cos(theta)
	cosA =  cA*np.cos(ksi) + sA*np.sin(ksi)
	sinA = -cA*np.sin(ksi) + sA*np.cos(ksi)
	alpha_det = np.arctan2(sinA,cosA)
	beta_det = np.arccos(cosB)
	return alpha_det, beta_det

#==================================//==================================//==============

def Normal_Vec(theta,phi): # Detector position in Geocentric Coordinates
	Zx = np.sin(theta)*np.cos(phi)
	Zy = np.sin(theta)*np.sin(phi)
	Zz = np.cos(theta)
	return Zx, Zy, Zz

def Nvec(alpha,beta,theta,phi,ksi): # Sky position in Detector Coordinates
	alpha_det, beta_det = AngTransf(alpha,beta,theta,phi,ksi)
	Nx = np.sin(beta_det)*np.cos(alpha_det)
	Ny = np.sin(beta_det)*np.sin(alpha_det)
	Nz = np.cos(beta_det)
	return Nx, Ny, Nz 

def Lvec0(alpha,beta,iota,psi): # Compute L in Geocentric Coordinates
	Lx = -np.sin(iota)*np.sin(psi)*np.cos(beta)*np.cos(alpha) - np.sin(iota)*np.cos(psi)*np.sin(alpha) + np.cos(iota)*np.sin(beta)*np.cos(alpha)
	Ly = -np.sin(iota)*np.sin(psi)*np.cos(beta)*np.sin(alpha) + np.sin(iota)*np.cos(psi)*np.cos(alpha) + np.cos(iota)*np.sin(beta)*np.sin(alpha)
	Lz = np.sin(iota)*np.sin(psi)*np.sin(beta) + np.cos(iota)*np.cos(beta)
	return Lx, Ly, Lz

def Lvec(L0,theta,phi,ksi): # Compute L in Detector Coordinates
	Lx0, Ly0, Lz0 = L0
	r11 = np.cos(theta)*np.cos(phi)*np.cos(ksi)-np.sin(phi)*np.sin(ksi)
	r12 = np.cos(theta)*np.sin(phi)*np.cos(ksi) + np.cos(phi)*np.sin(ksi)
	r13 = -np.sin(theta)*np.cos(ksi)
	r21 = -np.cos(theta)*np.cos(phi)*np.sin(ksi) - np.sin(phi)*np.cos(ksi)
	r22 = -np.cos(theta)*np.sin(phi)*np.sin(ksi) + np.cos(phi)*np.cos(ksi)
	r23 = np.sin(theta)*np.sin(ksi)
	r31 = np.sin(theta)*np.cos(phi)
	r32 = np.sin(theta)*np.sin(phi)
	r33 = np.cos(theta)
	Lxi = r11*Lx0 + r12*Ly0 + r13*Lz0
	Lyi = r21*Lx0 + r22*Ly0 + r23*Lz0
	Lzi = r31*Lx0 + r32*Ly0 + r33*Lz0
	return Lxi, Lyi, Lzi

def poll_ang(alpha0,beta0,iota,psi0,theta,phi,ksi):
	N = Normal_Vec(beta0,alpha0)
	Z = Normal_Vec(theta,phi)
	L = Lvec0(alpha0,beta0,iota,psi0)

	sPsi = scl_prod(Z,L) - scl_prod(N,L)*scl_prod(N,Z)
	cPsi = scl_prod(Z,vec_prod(N,L))
	psi = np.arctan2(sPsi,cPsi)
	return psi

def ObsAngles(alpha0,beta0,iota,psi0,lon,lat,rot):
	phi   = lon*np.pi/180
	theta = (90-lat)*np.pi/180
	ksi   = rot*np.pi/180

	alpha_obs, beta_obs = AngTransf(alpha0,beta0,theta,phi,ksi)
	psi_obs = poll_ang(alpha0,beta0,iota,psi0,theta,phi,ksi)
	return alpha_obs, beta_obs, psi_obs

rad = np.pi/180
R_earth = 6371.e3
c = 299792458.

def get_TimeDelay(det,prms):
	alpha0 = prms['RA']*rad
	beta0  = (90-prms['Dec'])*rad
	phi    = det['lon']*rad
	theta  = (90-det['lat'])*rad
	ksi    = det['rot']*rad

	_, beta_a = AngTransf(alpha0,beta0,theta,phi,ksi)
	return np.cos(beta_a)*R_earth/c