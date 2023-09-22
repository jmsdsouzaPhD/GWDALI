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

def pmt(vec):
	x, y = [ list(xx) for xx in list(permutations(vec)) ] , []
	for a in x:
		if(not a in y): y.append(a)
	return y

def Get_Tensors(FreeParams, gw_prms, detectors, dali_method, approximant, hide_info):
	
	Np = len(FreeParams)
	rho2 = 0

	Fisher   = np.zeros([Np, Np])
	Doublet3 = np.zeros([Np, Np, Np])
	Doublet4 = np.zeros([Np, Np, Np, Np])
	Triplet4 = np.zeros([Np, Np, Np, Np])
	Triplet5 = np.zeros([Np, Np, Np, Np, Np])
	Triplet6 = np.zeros([Np, Np, Np, Np, Np, Np])

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
		# Computing Doublet/Triplet (arXiv:2203.02670)
		#------------------------------------------------------
		
		idxs_doub, idxs_trip = sym.get_indep_indexies(Np)

		if(dali_method in ["Doublet","Triplet"]):
			indepD_12, indepD_22 = idxs_doub
			# --------------Doublet(1,2)--------------		
			for Idx2 in indepD_12:
				i, j, k = Idx2
				value = gwfunc.func_doublet3(*[FreeParams[ii] for ii in Idx2],gw_prms,det,approximant)
				for p in pmt([j,k]):
					jj, kk = p
					Doublet3[i][jj][kk] += value
			# --------------Doublet(2,2)--------------
			for Idx2 in indepD_22:
				i, j, k, l = Idx2
				value = gwfunc.func_doublet4(*[FreeParams[ii] for ii in Idx2],gw_prms,det,approximant)
				for p1 in pmt([i,j]):
					ii, jj = p1
					for p2 in pmt([k,l]):
						kk, ll = p2
						Doublet4[ii][jj][kk][ll] += value
			#*****************************#*****************************#*****************************
			if(dali_method == 'Triplet'):
				indepT_13, indepT_23, indepT_33 = idxs_trip
				# --------------Triplet(1,3)--------------
				for Idx3 in indepT_13:
					i, j, k, l = Idx3
					value = gwfunc.func_triplet4(*[FreeParams[ii] for ii in Idx3],gw_prms,det,approximant)
					for p in pmt([j,k,l]):
						jj, kk, ll = p
						Triplet4[i][jj][kk][ll] += value
				# --------------Triplet(2,3)--------------
				for Idx3 in indepT_23:
					i, j, k, l, m = Idx3
					value = gwfunc.func_triplet5(*[FreeParams[ii] for ii in Idx3],gw_prms,det,approximant)
					for p1 in pmt([i,j]):
						ii, jj = p1
						for p2 in pmt([k,l,m]):
							kk, ll, mm = p2
							Triplet5[ii][jj][kk][ll][mm] += value
				# --------------Triplet(3,3)--------------
				for Idx3 in indepT_33:
					i, j, k, l, m, n = Idx3
					value = gwfunc.func_triplet6(*[FreeParams[ii] for ii in Idx3],gw_prms,det,approximant)
					for p1 in pmt([i,j,k]):
						ii, jj, kk = p1
						for p2 in pmt([l,m,n]):
							ll, mm, nn = p2
							Triplet6[ii][jj][kk][ll][mm][nn] += value
							print('(%s) [Sym] Triplet: '%det['name'], i, j, k, l, m, n,'...'*10, end='\r')
		# --------------# --------------# --------------					

	SNR_sum = np.sqrt(rho2)
	if(not hide_info):
		print("\n >> SNR = ", SNR_sum)
		print(" >> 100/SNR = %.2e" % (100/SNR_sum) + '%\n')

	Doublet = [Doublet3.ravel(), Doublet4.ravel()]
	Triplet = [Triplet4.ravel(), Triplet5.ravel(), Triplet6.ravel()]

	return [SNR, GwData, Fisher , Doublet, Triplet]