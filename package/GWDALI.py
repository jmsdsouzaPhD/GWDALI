import numpy as np
import matplotlib.pyplot as plt
import os, shutil, sys
import warnings ; warnings.filterwarnings("ignore")
#===========================
#----Internal_libs----------
#===========================
import GWDALI.lib.Waveforms as wf
import GWDALI.lib.Angles_lib as geo
import GWDALI.lib.Dictionaries as gwdict
import GWDALI.lib.Diff_Signal_Tensors as gwfunc
import GWDALI.lib.Corner_Plots as corner
import GWDALI.lib.Auxiliar as Aux
from .lib.Likelihood import get_posterior
from .lib.Get_Tensors import Get_Tensors
#===========================
from scipy.interpolate import interp1d
from time import time as now
from tqdm import trange

#-------------------------------------------------
# From Dictionaries.py
PSD, labels_tex = gwdict.Load_Dictionaries()

#-------------------------------------------------------------------------------
#	"Compute the (Moore-Penrose) pseudo-inverse of a matrix"
#-------------------------------------------------------------------------------
def invert_matrix(matrix, rcond):

    dm = np.sqrt(np.diag(matrix))
    norm = np.outer(dm, dm)
    M_norm = matrix / norm

    try: M_inv = np.linalg.pinv(M_norm,rcond=rcond) / norm
    except: M_inv = np.linalg.pinv(matrix,rcond=rcond)

    return M_inv
#-------------------------------------------------------
# Get Waveforms and SNR
#-------------------------------------------------------
def get_Sn(det):
	name = det['name']
	try:
		Sn0 = PSD[name]
		freq0 = 10**np.linspace(0,4,4000)
		func_Sn = interp1d(freq0, Sn0, fill_value=np.inf, bounds_error=False)
	except:
		try:
			Sn0 = det['Sn']
			freq0 = det['freq']
		except:
			print("Fail on loading detectors!")
			quit()
	return Sn0, freq0

def get_strain_snr(gw_prms,detectors,approximant='TaylorF2',fmin=1.,fmax=1.e4,fsize=3.e3):
	Strains, SNR = [], []
	rho2 = 0

	freq = 10**np.linspace(np.log10(fmin), np.log10(fmax), int(fsize))
	for idx in range(len(detectors)):
		det = detectors[idx]
		Sn0, freq0 = get_Sn(det)
		func_Sn = interp1d(freq0,Sn0,fill_value=np.inf,bounds_error=False)
		detectors[idx]['freq'] = freq
		detectors[idx]['Sn']   = func_Sn(freq)

		h 	  = gwfunc.Signal(gw_prms, det, approximant)
		SNR2  = gwfunc.ScalarProduct(det['freq'], det['Sn'],h,h)
		rho2 += SNR2
		#print('#####'*10)
		Strains.append(h)
		SNR.append(np.sqrt(SNR2))

	hp, hx, freq = gwfunc.GW_Polarizations(gw_prms,freq, approximant)

	return [hp,hx,freq], Strains, SNR, np.sqrt(rho2)
#-------------------------------------------------------
# Main Function
#-------------------------------------------------------
def GWDALI( Detection_Dict, 
			FreeParams, 
			detectors, 
			approximant = 'TaylorF2_py',
			fmin  = 1., 
			fmax  = 1.e4, 
			fsize = 3000, 
			dali_method    = 'Doublet',
			sampler_method = 'nestle', 
			npoints      = 100,
			rcond		 = 1.e-4,
			new_priors   = None,
			save_samples = False, 
			save_cov     = False, 
			save_fisher  = False,
			plot_corner  = False,
			hide_info    = False,
			step_size	 = 1.e-6,
			diff_order	 = 2,
			run_sampler  = True,
			index		 = 1):

	path = Aux.CheckPath(save_cov, save_samples, save_fisher, plot_corner,dali_method)
	if(not hide_info):
		print('\n'+"-----"*10)
		print("-------------- Running GWDALI --------------------")
		print("-----"*10+'\n')
		print(" >> FreeParams: " ,  FreeParams)
		print(" >> dali_method: %s" % dali_method)
		print(" >> approximant: %s" % approximant)
		print(" >> sampler_method: %s\n" % sampler_method)

	Np   = len(FreeParams)
	keys = Detection_Dict.keys()
	freq = 10**np.linspace(np.log10(fmin), np.log10(fmax), fsize)

	#-------------------------------------------------
	# Setting Sn(f) by interpolation
	for idx in range(len(detectors)):
		det = detectors[idx]
		Sn0, freq0 = get_Sn(det)
		func_Sn = interp1d(freq0,Sn0,fill_value=np.inf,bounds_error=False)
		detectors[idx]['freq'] = freq
		detectors[idx]['Sn']   = func_Sn(freq)
	#-------------------------------------------------------
	# Store the injection values (Dict)
	truths = {}
	for key in Detection_Dict.keys():
		truths[key] = Detection_Dict[key]
	#-------------------------------------------------------
	# Store the injection values (List)
	Theta0 = np.array( [Detection_Dict[fp] for fp in FreeParams] )
	#-----------------------------------------------------------------------------------
	# From Compute_Tensors.py (looping in detectors)
	#-----------------------------------------------------------------------------------
	t1 = now() # initial time
	SNR , GwData, Fisher , Doublet , Triplet = Get_Tensors(FreeParams, Detection_Dict, detectors, dali_method, approximant, step_size, diff_order, hide_info)
	Tensors = [Fisher , Doublet , Triplet]
	t2 = now() # final time
	dT = t2 - t1
	dT_min = dT // 60
	dT_sec = dT % 60
	#-----------------------------------------------------------------------------------
	if(not hide_info): print("\n>> Fisher/Doublet/Triplet computed in %d sec" %dT)
	#------------------------------------------------------
	# Invert Fisher Matrix
	#------------------------------------------------------
	try: Cov = invert_matrix(Fisher, rcond) # from numpy.linalg.pinv
	except: Cov = np.ones([Np,Np])*np.nan
	
	Result = {}
	Result['Fisher']    = Fisher
	Result['CovFisher'] = Cov

	if(not hide_info):
		print("\n Fisher = \n", np.matrix(Fisher),'\n')
		print("\n Covariance = \n", np.matrix(Cov),'\n')
		print("-----"*10+'\n')

	if(save_fisher): Aux.Save_fisher(path,index,Detection_Dict,FreeParams,Fisher)
	if(save_cov): Aux.Save_CovFish(path,index,Detection_Dict,FreeParams,Cov)	
	#------------------------------------------------------
	# Monte-Carlo From Likelihood.py
	#------------------------------------------------------
	vec = [ truths[fp] for fp in FreeParams]
	truths = vec
	Cov2, M, Error, Recovery = [],[],[],[]
	if(dali_method != 'Fisher'):
		# allowed only for Fisher_Sampling and Doublet
		if(run_sampler):
			if(not hide_info): print("Running Sampler ...\n")
			M = get_posterior(FreeParams, Theta0, Detection_Dict, GwData, approximant, detectors, Tensors, dali_method, sampler_method, npoints, new_priors)
			qs = [0.2,0.5,0.8] # CL = 60%
			for m in M.T:
				q = np.quantile(m,qs)
				Recovery.append(q[1])
				Error.append(0.5*(q[2]-q[0]))
			Cov2 = np.matrix( np.cov(np.transpose(M)) )
			if(save_cov): Aux.Save_Cov(path,index,Detection_Dict,FreeParams,M,Cov2)
				
		Result['Covariance'] = Cov2
		Result['Samples']    = M
		Result['Error']      = Error
		Result['Recovery']   = Recovery
	
	Result['SNR']	  = SNR
	Result['Tensors'] = Tensors

	if(save_samples): Aux.Save_Samples(path,index,Detection_Dict,FreeParams,M)
	if(plot_corner): Aux.PlotCorner(M,truths,Cov,FreeParams,[dT_min,dT_sec],dali_method,path,index)

	# Result (dict)
	# keys: Fisher, CovFisher, Covariance, Samples, Error, Recovery, SNR
	return Result 