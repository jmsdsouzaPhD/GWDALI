
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
import GWDALI.lib.Derivatives_Tensors as gwfunc
import GWDALI.lib.Corner_Plots as corner

from .lib.Likelihood import get_posterior
from .lib.Compute_Tensors import Get_Tensors

#===========================

from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
from time import time as now
from tqdm import trange

cosmo = FlatLambdaCDM(70,0.3)

#-------------------------------------------------
# From Dictionaries.py
PSD, labels_tex = gwdict.Load_Dictionaries()
#-------------------------------------------------
freq0 = 10**np.linspace(0,4,4000)

#-------------------------------------------------------------------------------
#	"Compute the (Moore-Penrose) pseudo-inverse of a matrix"
#-------------------------------------------------------------------------------
def invert_matrix(matrix, rcond): 

    dm = np.sqrt(np.diag(matrix))
    norm = np.outer(dm, dm)
    M_norm = matrix / norm

    M_inv = np.linalg.pinv(M_norm,rcond=rcond)

    return M_inv / norm
#-------------------------------------------------------
# Get Waveforms and SNR
#-------------------------------------------------------
def get_strain_snr(gw_prms,detectors,approximant='TaylorF2',fmin=1.,fmax=1.e4,fsize=3.e3):
	Strains, SNR = [], []
	rho2 = 0

	freq = 10**np.linspace(np.log10(fmin), np.log10(fmax), int(fsize))
	for idx in range(len(detectors)):
		det = detectors[idx]
		name = det['name']
		Sn0 = PSD[name]
		func_Sn = interp1d(freq0, Sn0, fill_value=np.inf, bounds_error=False)
		detectors[idx]['freq'] = freq
		detectors[idx]['Sn']   = func_Sn(freq)

		h 	  = gwfunc.Signal(gw_prms, det, approximant)
		SNR2  = gwfunc.ScalarProduct(det['freq'], det['Sn'],h,h)
		rho2 += SNR2
		print('#####'*10)
		Strains.append(h)
		SNR.append(np.sqrt(SNR2))

	hp, hx = gwfunc.GW_Polarizations(gw_prms,freq, approximant)

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
			index		 = 1):

	if(not hide_info):
		print('\n'+"-----"*10)
		print("-------------- Running GWDALI --------------------")
		print("-----"*10+'\n')
		print(" >> FreeParams: " ,  FreeParams)
		print(" >> dali_method: %s" % dali_method)
		print(" >> approximant: %s" % approximant)
		print(" >> sampler_method: %s\n" % sampler_method)

	T1 = now()
	
	cond_out = any([save_cov, save_samples, save_fisher, plot_corner])
	#-----------------------------------------------------------------------------------
	# Make a directory to store the results
	if(cond_out):
		path = 'output_%s/' % dali_method
		if(not os.path.isdir(path)): os.mkdir(path)
	#-----------------------------------------------------------------------------------

	Np   = len(FreeParams)
	keys = Detection_Dict.keys()
	freq = 10**np.linspace(np.log10(fmin), np.log10(fmax), fsize)

	#-------------------------------------------------
	# Setting Sn(f) by interpolation
	for idx in range(len(detectors)):
		det = detectors[idx]
		name = det['name']
		Sn0 = PSD[name]
		func_Sn = interp1d(freq0, Sn0, fill_value=np.inf, bounds_error=False)
		detectors[idx]['freq'] = freq
		detectors[idx]['Sn']   = func_Sn(freq)

	t_init = now()
	#-----------------------------------------------------------------------------------
	
	Result = {}

	#-------------------------------------------------------
	# Store the injection values (Dict)
	truths = {}
	for key in Detection_Dict.keys():
		truths[key] = Detection_Dict[key]
	#-------------------------------------------------------
	# Store the injection values (List)
	Theta0 = []
	for fp in FreeParams:
		Theta0.append(Detection_Dict[fp])
	Theta0 = np.array(Theta0)
	#-------------------------------------------------------
	
	t1 = now()
	#-----------------------------------------------------------------------------------
	# From Compute_Tensors.py (looping in detectors)
	#-----------------------------------------------------------------------------------
	SNR , GwData, Fisher , Doublet , Triplet = Get_Tensors(FreeParams, Detection_Dict, detectors, dali_method, approximant,hide_info)
	Tensors = [Fisher , Doublet , Triplet]
	#-----------------------------------------------------------------------------------
	t2 = now()
	dt = t2 - t1
	if(not hide_info):
		print('\n'+"-----"*10)
		print("\n>> Fisher/Doublet/Triplet computed in %d sec" %dt)
		print('\n'+"-----"*10)

	#********************************************************************

	#------------------------------------------------------
	# Invert Fisher Matrix
	#------------------------------------------------------
	
	Cov = invert_matrix(Fisher, rcond) # from numpy.linalg.pinv

	Result['Fisher']    = Fisher
	Result['CovFisher'] = Cov

	if(not hide_info):
		print("\n Fisher = \n", np.matrix(Fisher),'\n')
		print("\n Covariance = \n", np.matrix(Cov),'\n')
		print("-----"*10+'\n')

	#------------------------------------------------------
	# Save the injection in each file headers:
	#------------------------------------------------------
	header01 = 'FullInjections:\n'
	header02 = ''
	for key in keys:
		header01 += "%s\t" % key
		header02 += "%e\t" % Detection_Dict[key]
	header1 = 'Injections: ' ; header2 = '' ; line = '---'*10
	for fp in FreeParams:
		header1 += '%s\t' % fp
		header2 += '%e\t' % Detection_Dict[fp]
	header = header01 + '\n' + header02 + '\n' + line + '\n' + header1+'\n'+header2 +'\n'+ line

	if(save_fisher): np.savetxt(path+'Fisher_Matrix_%d.txt'%(index),Fisher,fmt='%e',delimiter='\t',header=header)
	if(save_cov): np.savetxt(path+'FisherCov_%d.txt'%(index),Cov,fmt='%e',delimiter='\t',header=header)

	#------------------------------------------------------#------------------------------------------------------
	#-------------------------------------------------MONTE-CARLO-------------------------------------------------
	#------------------------------------------------------#------------------------------------------------------
	vec = [ truths[fp] for fp in FreeParams]
	truths = vec		
	#------------------------------------------------------
	# From Likelihood.py
	#------------------------------------------------------
	if(dali_method != 'Fisher'):
		# allowed only for Fisher_Sampling and Doublet
		if(not hide_info): print("Running Sampler ...\n")
		#------------------------------------------------------#------------------------------------------------------
		M = get_posterior(FreeParams, Theta0, Detection_Dict, GwData, approximant, detectors, Tensors, dali_method, sampler_method, npoints, new_priors)
		#------------------------------------------------------#------------------------------------------------------
		Cov2 = np.matrix( np.cov(np.transpose(M)) )
		header3 = header + '\n' + line + '\nRecovery:\n'
		for i in range(len(FreeParams)): header3 += '%e\t' % np.average(M[:,i])
		header3 += '\n' + line
		if(save_cov): np.savetxt(path+'Covariance_Matrix_%d.txt'%(index),Cov2,fmt='%e',delimiter='\t',header=header3)

		qs = [0.2,0.5,0.8] # CL = 60%
		Error, Recovery = [], []
		for m in M.T:
			q = np.quantile(m,qs)
			Recovery.append(q[1])
			Error.append(0.5*(q[2]-q[0]))

		Result['Covariance'] = Cov2
		Result['Samples']    = M
		Result['Error']      = Error
		Result['Recovery']   = Recovery
	Result['SNR']		 = SNR
	#------------------------------------------------------
	# Save Samples
	#------------------------------------------------------

	T2 = now()
	dT = T2 - T1
	dT_min = dT // 60
	dT_sec = dT % 60
	title = "Concluded in [%d min : %d sec]" % (dT_min, dT_sec)

	if(save_samples and (dali_method!='Fisher')):
		if(not hide_info): print(">> Saving Samples ...")
		np.savetxt(path+'samples_%d.txt' % (index),M,fmt='%e',delimiter='\t',header=header)

	#------------------------------------------------------
	# Plot Corner
	#------------------------------------------------------

	if(plot_corner):
		labels = [ labels_tex[fp] for fp in FreeParams]

		if(dali_method == 'Fisher'):
			fig = corner.Plot_Fisher(truths, labels, FreeParams, Cov)
			plt.suptitle(title)
			fig.subplots_adjust(wspace=0.01, hspace=0.01)
			fig.savefig(path+'corner_%d.png' % (index))
		else:
			print(">>> Plotting Corner [Doublet] .... ")
			fig = corner.Plot_corner(M, truths, labels, FreeParams, Cov)
			plt.suptitle(title)
			fig.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95))
			fig.subplots_adjust(wspace=0.01, hspace=0.01)
			fig.savefig(path+'corner_%d.png' % (index))	
		plt.show()

	# Result (dict)
	# keys: Fisher, CovFisher, Covariance, Samples, Error, Recovery, SNR
	return Result 