import numpy as np
from scipy.interpolate import interp1d, CubicSpline

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

M_sun = 1.98e30
G = 6.673e-11
c = 299792458. # m/s

#==============================================================================================
# LAL WAVEFORMS
#==============================================================================================

def Waveform_lal(m1, m2, iota, DL, s1, s2, freq, approx):
	dF = np.average(np.diff(freq)) #1 
	f_low, f_max = min(freq), max(freq)

	sx1, sy1, sz1 = s1
	sx2, xy2, sz2 = s2
	
	phi_ref = 0.
	f_ref = 0.

	DL *= 3.086e22 # Mpc -> meters
	hp,hx=lalsim.SimInspiralChooseFDWaveform(m1*M_sun,m2*M_sun, sx1, sy1, sz1, sx2, xy2, sz2, DL,
											iota, phi_ref, 0., 0., 0., dF, f_low, f_max, f_ref, None,
											lalsim.GetApproximantFromString(approx))
	
	hp = hp.data.data ; hp[np.isnan(hp)] = 0
	hx = hx.data.data ; hx[np.isnan(hx)] = 0

	N = len(hp)
	freq0 = np.linspace(0,dF*(N-1),N)

	return hp, hx, freq0

#==============================================================================================
# PYTHON WAVEFORMS ( Leading Order and TaylorF2(3.5PN) )
#==============================================================================================
	
def GW_Phase_TaylorF2(m1,m2,freq):
	M = m1+m2
	eta = m1*m2/M**2
	M *= M_sun
	GM_c3 = G*M/c**3
	
	A = list(np.zeros(8))
	A[0] = 1.
	A[1] = 0.
	A[2] = (3715./756) + (55./9)*eta 
	A[3] = -16*np.pi
	A[4] = (15293365./508032) + (27145*eta/504) + (3085*eta**2/72)
	A[5] = np.pi*((38645./756) - (65*eta/9))*(1.+np.log((6.**1.5)*np.pi*GM_c3*freq))
	A[6] = (11583231236531./4694215680) - (640*np.pi**2/3) - (6848*np.euler_gamma/21)\
		 + (-(15737765635./3048192)+(2255*np.pi**2/12))*eta + (76055*eta**2/1728)\
		 - (127825*eta**3/1296) - (6848./63)*np.log(64*np.pi*GM_c3*freq)
	A[7] = np.pi*((77096675./254016) + (378515*eta/1512) - (74045*eta**2/756))

	phase = np.zeros(len(freq))  # 3.5PN Corrections from idx>0
	for idx in range(0,8):
		phase += A[idx]*(np.pi*GM_c3*freq)**((idx-5)/3)
	phase *= 3./(128*eta)

	return phase

def GW_Amplitude(m1,m2,DL,freq):
	DL *= 3.086e22     # DL [meters]
	M = m1+m2			 # Total Mass
	eta = m1*m2/M**2   # Symmetric Mass Ratio
	M *= M_sun
	Mc = M*eta**(3./5) # Chirp Mass
	R_isco = 6*(G*M)/c**2
	f_isco = np.sqrt(G*M/R_isco**3)/np.pi # gw frequency at R=R_isco
	Cutoff = np.ones(len(freq))*(freq<4*f_isco)
	return Cutoff*pow(freq,-7./6)*np.sqrt(5./(24*c**3))*pow(np.pi,-2./3)*pow(G*Mc,5./6)/DL

def Waveform_TaylorF2(m1, m2, iota, DL, s1, s2, freq): # Sathyaprakash-Schutz(2009)
	Amp   = GW_Amplitude(m1,m2,DL,freq)
	phase = GW_Phase_TaylorF2(m1,m2,freq)
	
	exp = np.exp(1.j*(phase-np.pi/4))

	u = np.cos(iota)
	Gp, Gx = 0.5*(1.+u*u), 1.j*u
	Hplus  = Gp*Amp*exp
	Hcross = Gx*Amp*exp

	return Hplus, Hcross, freq
		
def Waveform_Simple(m1, m2, iota, DL, s1, s2, freq):
	m1 *= M_sun ; m2 *= M_sun
	M = m1+m2		   # Total Mass
	eta = m1*m2/M**2   # Symmetric Mass Ratio
	Mc = eta**(3./5) * M
	DL *= 3.086e22

	h0 = np.sqrt(5./24) * (G*Mc/c**3)**(5./6) / (np.pi**(2./3)*DL/c)
	phase = (3./128) * (np.pi * G * Mc * freq / c**3)**(-5./3)

	R_isco = 6*G*M/c**2
	f_isco = np.sqrt(G*M/R_isco**3) / np.pi
	cutoff = freq<f_isco

	Phase = - np.pi/4 + (3./128)*(np.pi*G*Mc*freq/c**3)**(-5./3)
	u = np.cos(iota)

	Gp = 0.5*(1+u*u)
	Gx = 1.j*u

	hp = h0*Gp*freq**(-7./6)*np.exp(1.j*Phase) * cutoff
	hx = h0*Gx*freq**(-7./6)*np.exp(1.j*Phase) * cutoff
	
	return hp, hx, freq

#-----------------------------------------------------------------------------------------------

def Waveforms(m1, m2, iota, DL, s1, s2, freq, approx):
	if approx == 'Leading_Order':
		return Waveform_Simple(m1, m2, iota, DL, s1, s2, freq)
	elif approx == 'TaylorF2_py':
		return Waveform_TaylorF2(m1, m2, iota, DL, s1, s2, freq)
	else:
		return Waveform_lal(m1, m2, iota, DL, s1, s2, freq, approx)