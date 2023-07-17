import os, sys
import numpy as np
import GWDALI.lib.Angles_lib as geo
import GWDALI.lib.Dictionaries as gwdict
import GWDALI.lib.Derivatives_Tensors as gwfunc
import GWDALI.lib.Symmetric_Indexies as sym

from itertools import permutations
from tqdm import trange

deg2rad = np.pi/180

PSD, labels_tex = gwdict.Load_Dictionaries()

def get_idx(idxs,Np):
	n = len(idxs)-1 ; IdxP = 0
	for i in range(len(idxs)): 
		IdxP += idxs[i]*Np**n ; n-=1
	return IdxP

def pmt(vec):
	x, y = [ list(xx) for xx in list(permutations(vec)) ] , []
	for a in x:
		if(not a in y): y.append(a)
	return y

def Get_Tensors(FreeParams, gw_prms, detectors, dali_method, approximant, hide_info):
	
	#if(hide_info): sys.stdout = open(os.devnull, 'w')

	Np = len(FreeParams)
	rho2 = 0

	Fisher   = np.zeros([Np, Np])
	Doublet3 = np.zeros(Np**3)
	Doublet4 = np.zeros(Np**4)
	Triplet4 = np.zeros(Np**4)
	Triplet5 = np.zeros(Np**5)
	Triplet6 = np.zeros(Np**6)

	num_det = 1
	GwData = [] ; SNR = []
	for det in detectors:

		#-----------------------------------
		# Computing Signal to Noise 
		#-----------------------------------
		# from Derivative_Tensors.py
		h 	  = gwfunc.Signal(gw_prms, det, approximant)
		SNR2  = gwfunc.ScalarProduct(det['freq'], det['Sn'],h,h)
		rho2 += SNR2
		GwData.append(h) ; SNR.append(np.sqrt(SNR2))
		#------------------------------------------------------
		# Computing Fisher Matrix
		#------------------------------------------------------

		if(not hide_info): print("\n\n\t Computing Fisher %s(%d) ..." % (det['name'],num_det))
		for i in range(Np):
			for j in range(i+1):
				if(not hide_info): print('Fisher: [%d-%d]' % (i+1,j+1), end='\r')
				xj = FreeParams[j]
				xi = FreeParams[i]

				Fij = gwfunc.Fisher_ij(xi, xj, gw_prms, det, approximant)
				Fisher[i][j] += Fij
				if(i!=j): Fisher[j][i] += Fij

		#------------------------------------------------------
		# Computing Doublet/Triplet
		#------------------------------------------------------
		
		idxs_doub, idxs_trip = sym.get_indep_indexies(Np)

		if(dali_method in ["Doublet","Triplet"]):
			indepD_12, indepD_22 = idxs_doub
			# --------------Doublet(1,2)--------------		
			for Idx2 in indepD_12:
				i, j, k = Idx2
				value = gwfunc.func_doublet3(*[FreeParams[ii] for ii in Idx2],gw_prms,det,approximant)
				for p in pmt([j,k]):
					Doublet3[get_idx([i]+p,Np)] += value
			# --------------Doublet(2,2)--------------
			for Idx2 in indepD_22:
				i, j, k, l = Idx2
				value = gwfunc.func_doublet4(*[FreeParams[ii] for ii in Idx2],gw_prms,det,approximant)
				for p1 in pmt([i,j]):
					for p2 in pmt([k,l]):
						Doublet4[get_idx(p1+p2,Np)] += value
			#*****************************#*****************************#*****************************
			if(dali_method == 'Triplet'):
				indepT_13, indepT_23, indepT_33 = idxs_trip
				# --------------Triplet(1,3)--------------
				for Idx3 in indepT_13:
					i, j, k, l = Idx3
					value = gwfunc.func_triplet4(*[FreeParams[ii] for ii in Idx3],gw_prms,det,approximant)
					for p in pmt([j,k,l]):
						Triplet4[get_idx([i]+p,Np)] += value
				# --------------Triplet(2,3)--------------
				for Idx3 in indepT_23:
					i, j, k, l, m = Idx3
					value = gwfunc.func_triplet5(*[FreeParams[ii] for ii in Idx3],gw_prms,det,approximant)
					for p1 in pmt([i,j]):
						for p2 in pmt([k,l,m]):
							Triplet5[get_idx(p1+p2,Np)] += value
				# --------------Triplet(3,3)--------------
				for Idx3 in indepT_33:
					i, j, k, l, m, n = Idx3
					value = gwfunc.func_triplet6(*[FreeParams[ii] for ii in Idx3],gw_prms,det,approximant)
					for p1 in pmt([i,j,k]):
						for p2 in pmt([l,m,n]):
							Triplet6[get_idx(p1+p2,Np)] += value
							print('(%s) [Sym] Triplet: '%det['name'], i, j, k, l, m, n,'...'*10, end='\r')
		# --------------# --------------# --------------					
	Doublet = [Doublet3, Doublet4]
	Triplet = [Triplet4, Triplet5, Triplet6]
	SNR_sum = np.sqrt(rho2)
	if(not hide_info):
		print("\n >> SNR = ", SNR_sum)
		print(" >> 100/SNR = %.2e" % (100/SNR_sum) + '%\n')

	Doublet = [Doublet3, Doublet4]
	Triplet = [Triplet4, Triplet5, Triplet6]

	return [SNR, GwData, Fisher , Doublet, Triplet]