import os, sys
import numpy as np
import GWDALI.lib.Angles_lib as geo
import GWDALI.lib.Dictionaries as dct
import GWDALI.lib.Diff_Signal_Tensors as gw
import GWDALI.lib.Auxiliar as sym

from itertools import permutations
from tqdm import trange

deg2rad = np.pi/180

PSD, labels_tex = dct.Load_Dictionaries()

def pmt(vec):
	x, y = [ list(xx) for xx in list(permutations(vec)) ] , []
	for a in x:
		if(not a in y): y.append(a)
	return y

def Get_Tensors(FreeParams, gw_prms, detectors, dali_method, approx, step_size, diff_order, hide_info):
	
	Np = len(FreeParams)
	rho2 = 0

	Fisher = np.zeros([Np, Np])
	Db3    = np.zeros([Np, Np, Np])
	Db4    = np.zeros([Np, Np, Np, Np])
	Tp4    = np.zeros([Np, Np, Np, Np])
	Tp5    = np.zeros([Np, Np, Np, Np, Np])
	Tp6    = np.zeros([Np, Np, Np, Np, Np, Np])

	num_det = 1
	GwData = [] ; SNR = []
	for det in detectors:

		#-----------------------------------
		# Computing Signal to Noise 
		#-----------------------------------
		# from Derivative_Tensors.py
		h 	  = gw.Signal(gw_prms, det, approx)
		SNR2  = gw.ScalarProduct(det['freq'], det['Sn'],h,h)
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
				Fij = gw.Fisher_ij(xi, xj, gw_prms, det, approx, step_size, diff_order)
				Fisher[i][j] += Fij

				if(i!=j): Fisher[j][i] += Fij
		if(not hide_info): print("\n\n\t\tFisher Concluded!!!\n\n")
		#------------------------------------------------------
		# Computing Doublet/Triplet (arXiv:2203.02670)
		#------------------------------------------------------
		
		idxs_doub, idxs_trip = sym.get_indep_indexies(Np)
		
		if(dali_method in ["Doublet","Triplet"]):
			indepD_12, indepD_22 = idxs_doub
			# --------------Doublet(1,2)--------------		
			for Idx2 in indepD_12:
				i, j, k = Idx2
				value = gw.Doublet3(*[FreeParams[ii] for ii in Idx2],gw_prms,det,approx, step_size, diff_order)
				for p in pmt([j,k]):
					jj, kk = p
					Db3[i][jj][kk] += value
			# --------------Doublet(2,2)--------------
			for Idx2 in indepD_22:
				i, j, k, l = Idx2
				value = gw.Doublet4(*[FreeParams[ii] for ii in Idx2],gw_prms,det,approx, step_size, diff_order)
				for p1 in pmt([i,j]):
					ii, jj = p1
					for p2 in pmt([k,l]):
						kk, ll = p2
						Db4[ii][jj][kk][ll] += value
						if(not hide_info): print('(%s) [Sym] Doublet: '%det['name'], i, j, k, l,'...'*10, end='\r')
			#*****************************#*****************************#*****************************
			if(dali_method == 'Triplet'):
				indepT_13, indepT_23, indepT_33 = idxs_trip
				# --------------Triplet(1,3)--------------
				for Idx3 in indepT_13:
					i, j, k, l = Idx3
					value = gw.Triplet4(*[FreeParams[ii] for ii in Idx3],gw_prms,det,approx, step_size, diff_order)
					for p in pmt([j,k,l]):
						jj, kk, ll = p
						Tp4[i][jj][kk][ll] += value
				# --------------Triplet(2,3)--------------
				for Idx3 in indepT_23:
					i, j, k, l, m = Idx3
					value = gw.Triplet5(*[FreeParams[ii] for ii in Idx3],gw_prms,det,approx, step_size, diff_order)
					for p1 in pmt([i,j]):
						ii, jj = p1
						for p2 in pmt([k,l,m]):
							kk, ll, mm = p2
							Tp5[ii][jj][kk][ll][mm] += value
				# --------------Triplet(3,3)--------------
				for Idx3 in indepT_33:
					i, j, k, l, m, n = Idx3
					value = gw.Triplet6(*[FreeParams[ii] for ii in Idx3],gw_prms,det,approx, step_size, diff_order)
					for p1 in pmt([i,j,k]):
						ii, jj, kk = p1
						for p2 in pmt([l,m,n]):
							ll, mm, nn = p2
							Tp6[ii][jj][kk][ll][mm][nn] += value
							if(not hide_info): print('(%s) [Sym] Triplet: '%det['name'], i, j, k, l, m, n,'...'*10, end='\r')
		# --------------# --------------# --------------					

	SNR_sum = np.sqrt(rho2)
	if(not hide_info):
		print("\n >> SNR = ", SNR_sum)
		print(" >> 100/SNR = %.2e" % (100/SNR_sum) + '%\n')

	Doublet = [Db3, Db4]
	Triplet = [Tp4, Tp5, Tp6]

	return [SNR, GwData, Fisher , Doublet, Triplet]