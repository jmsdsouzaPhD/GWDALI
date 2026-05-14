import os
import numpy as np
import matplotlib.pyplot as plt
from itertools import permutations
import GWDALI.lib.Dictionaries as gwdict
import GWDALI.lib.Corner_Plots as corner

line = "---"*30
PSD, labels_tex = gwdict.Load_Dictionaries()

def Load_Header(Detection_Dict, FreeParams):
	hd1 = 'FullInjections:\n'
	hd2 = ''
	keys = Detection_Dict.keys()
	for key in keys:
		hd1 += "%s\t" % key
		hd2 += "%e\t" % Detection_Dict[key]
	header1 = 'Injections: ' ; header2 = '' ; line = '---'*10
	for fp in FreeParams:
		header1 += '%s\t' % fp
		header2 += '%e\t' % Detection_Dict[fp]
	header = hd1 + '\n' + hd2 + '\n' + line + '\n' + header1+'\n'+header2 +'\n'+ line
	return header

def CheckPath(save_cov, save_samples, save_fisher, plot_corner,dali_method):
	cond_out = any([save_cov, save_samples, save_fisher, plot_corner])
	path = None
	if(cond_out):
		path = 'output_%s/' % dali_method
		if(not os.path.isdir(path)): os.mkdir(path)
	return path

def Save_fisher(path,index,Detection_Dict,FreeParams,Fisher):
	header = Load_Header(Detection_Dict, FreeParams) 
	np.savetxt(path+'Fisher_Matrix_%d.txt'%(index),Fisher,fmt='%e',delimiter='\t',header=header)
	pass

def Save_CovFish(path,index,Detection_Dict,FreeParams,Cov):
	header = Load_Header(Detection_Dict, FreeParams)	
	np.savetxt(path+'FisherCov_%d.txt'%(index),Cov,fmt='%e',delimiter='\t',header=header)
	pass

def Save_Cov(path,index,Detection_Dict,FreeParams,M,Cov2):
	header = Load_Header(Detection_Dict, FreeParams)
	header3 = header + '\n' + line + '\nRecovery:\n'
	for i in range(len(FreeParams)): header3 += '%e\t' % np.average(M[:,i])
	header3 += '\n' + line
	np.savetxt(path+'Covariance_Matrix_%d.txt'%(index),Cov2,fmt='%e',delimiter='\t',header=header3)
	pass

def Save_Samples(path,index,Detection_Dict,FreeParams,M):
	header = Load_Header(Detection_Dict, FreeParams)
	np.savetxt(path+'samples_%d.txt' % (index),M,fmt='%e',delimiter='\t',header=header)

def PlotCorner(M,truths,Cov,FreeParams,dT,dali_method,path,index):
	dT_min, dT_sec = dT
	labels = [ labels_tex[fp] for fp in FreeParams]
	title = "Concluded in [%d min : %d sec]" % (dT_min, dT_sec)
	if(dali_method == 'Fisher'):
		fig = corner.Plot_Fisher(truths, labels, FreeParams, Cov)
		plt.suptitle(title)
		fig.subplots_adjust(wspace=0.01, hspace=0.01)
		fig.savefig(path+'corner_%d.png' % (index))
	else:
		fig = corner.Plot_corner(M, truths, labels, FreeParams, Cov)
		plt.suptitle(title)
		fig.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95))
		fig.subplots_adjust(wspace=0.01, hspace=0.01)
		fig.savefig(path+'corner_%d.png' % (index))	
	plt.show()
	pass

#===========================================#===========================================
#===========================================#===========================================

def get_indep_indexies(Np):
	#---------------------#---------------------
	# --------------Doublet(1,2)--------------
	#---------------------#---------------------
	indepD_12 = []
	for i in range(Np):
		for j in range(Np):
			for k in range(j,Np):
				indepD_12.append([i,j,k])
	#---------------------#---------------------
	# --------------Doublet(2,2)--------------
	#---------------------#---------------------
	indepD_22 = []
	for i in range(Np):
		for j in range(i,Np):
			for k in range(Np):
				for l in range(k,Np):
					indepD_22.append([i,j,k,l])
	#---------------------#---------------------
	# --------------Triplet(1,3)--------------
	#---------------------#---------------------
	indepT_13 = []
	for i in range(Np):
		for j in range(Np):
			for k in range(j,Np):
				for l in range(k,Np):
					indepT_13.append([i,j,k,l])
	#---------------------#---------------------
	# --------------Triplet(2,3)--------------
	#---------------------#---------------------
	indepT_23 = []
	for i in range(Np):
		for j in range(i,Np):
			for k in range(Np):
				for l in range(k,Np):
					for m in range(l,Np):
						indepT_23.append([i,j,k,l,m])
	#---------------------#---------------------
	# --------------Triplet(3,3)--------------
	#---------------------#---------------------
	indepT_33 = []
	for i in range(Np):
		for j in range(i,Np):
			for k in range(j,Np):
				for l in range(Np):
					for m in range(l,Np):
						for n in range(m,Np):
							indepT_33.append([i,j,k,l,m,n])
	#-----------------------------#-----------------------------
	idxs_doub = [indepD_12, indepD_22]
	idxs_trip = [indepT_13, indepT_23, indepT_33]

	return idxs_doub, idxs_trip