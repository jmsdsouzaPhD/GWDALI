import numpy as np
from scipy.interpolate import interp1d
from pathlib import Path

import GWDALI.lib.Waveforms as wf

#------------------------------------------------------
# PATTERN FUNCTIONS AND WAVEFORMS
#------------------------------------------------------

dets = ['ET','CE','aLIGO','aVirgo','KAGRA']

Waveforms = {}
Waveforms['Leading_Order'] = wf.Waveform_Simple
Waveforms['TaylorF2']      = wf.Waveform_TaylorF2
try:
	Waveforms['TaylorF2_lal']  = wf.Waveform_TaylorF2_lal
	Waveforms['IMRPhenomD']    = wf.Waveform_IMRPhenomD
	Waveforms['IMRPhenomP']    = wf.Waveform_IMRPhenomP
except:
	pass

#------------------------------------------------------
# DETECTORS SENSITIVITY
#------------------------------------------------------

path = Path(__file__).parent / '../Detectors_Sensitivity/'

data_LIGO  = np.loadtxt(path / 'Sn_L.txt')
data_Virgo = np.loadtxt(path / 'Sn_V.txt')
data_Kagra = np.loadtxt(path / 'Sn_K.txt')
data_ET    = np.loadtxt(path / 'Sn_ET.txt')
data_CE    = np.loadtxt(path / 'Sn_CE.txt')

freq1, Sn1	= data_LIGO[:,0]  , data_LIGO[:,1]**2
freq2, Sn2	= data_Virgo[:,0] , data_Virgo[:,1]**2
freq3, Sn3	= data_Kagra[:,0] , data_Kagra[:,1]**2
freq4, Sn4	= data_ET[:,0]	  , data_ET[:,1]**2
freq5, Sn5	= data_CE[:,0]	  , data_CE[:,1]**2

freq1, Sn1	= data_LIGO[:,0]  , data_LIGO[:,1]**2
freq2, Sn2	= data_Virgo[:,0] , data_Virgo[:,1]**2
freq3, Sn3	= data_Kagra[:,0] , data_Kagra[:,1]**2
freq4, Sn4	= data_ET[:,0]	  , data_ET[:,1]**2
freq5, Sn5	= data_CE[:,0]	  , data_CE[:,1]**2

func1 = interp1d(freq1, Sn1, fill_value=np.inf, bounds_error=False)
func2 = interp1d(freq2, Sn2, fill_value=np.inf, bounds_error=False)
func3 = interp1d(freq3, Sn3, fill_value=np.inf, bounds_error=False)
func4 = interp1d(freq4, Sn4, fill_value=np.inf, bounds_error=False)
func5 = interp1d(freq5, Sn5, fill_value=np.inf, bounds_error=False)

#------------------------------------------------------

freq0 = 10**np.linspace(0,4,4000)

SnL 	   = func1(freq0)
SnV, SnK   = func2(freq0), func3(freq0)
SnET, SnCE = func4(freq0), func5(freq0)

PSD = {
'aLIGO': SnL,
'aVirgo': SnV,
'KAGRA': SnK,
'CE': SnCE,
'ET': SnET
}

#------------------------------------------------------
# LATEX LABELS
#------------------------------------------------------

labels_tex = {} 
labels_tex['m1']       = '$m_1$ [$M_{\odot}$]'
labels_tex['m2']       = '$m_2$ [$M_{\odot}$]'
labels_tex['Mc']       = '$M_c$ [$M_{\odot}$]'
labels_tex['eta']      = '$\eta$'
labels_tex['q']        = '$q$'
labels_tex['DL']       = '$d_L$ [Gpc]'
labels_tex['iota']     = '$\iota$ [rad]'
labels_tex['psi']      = '$\psi$ [rad]'
labels_tex['alpha']    = '$\\alpha$ [rad]'
labels_tex['beta']     = '$\\beta$ [rad]'
labels_tex['t_coal']   = '$t_{coal}$ [sec]'
labels_tex['phi_coal'] = '$\phi_{coal}$ [rad]'
labels_tex['sx1']      = '$S_{x,1}$'
labels_tex['sy1']      = '$S_{y,1}$'
labels_tex['sz1']      = '$S_{z,1}$'
labels_tex['sx2']      = '$S_{x,2}$'
labels_tex['sy2']      = '$S_{y,2}$'
labels_tex['sz2']      = '$S_{z,2}$'

labels_tex['RA']  = '$RA$ [deg]'
labels_tex['Dec'] = '$Dec$ [deg]'

def Load_Dictionaries():
	return Waveforms, PSD, labels_tex

#------------------------------------------------------