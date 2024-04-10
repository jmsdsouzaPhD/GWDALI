import os, sys
import numpy as np
import GWDALI.lib.Dictionaries as dct
import GWDALI.lib.Diff_Signal_Tensors as gw
import GWDALI.lib.Auxiliar as sym

from itertools import permutations
from tqdm import trange

deg2rad = np.pi/180
c = 299792458. # m/s
R_earth = 6371.e3 # meters

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
	Db12    = np.zeros([Np, Np, Np])
	Db22    = np.zeros([Np, Np, Np, Np])
	Tp13    = np.zeros([Np, Np, Np, Np])
	Tp23    = np.zeros([Np, Np, Np, Np, Np])
	Tp33    = np.zeros([Np, Np, Np, Np, Np, Np])

	num_det = 1
	GwData = [] ; SNR = []
	det_a = detectors[0]
	for det_b in detectors:
		dets = [det_a, det_b]
		#-----------------------------------
		# Computing Signal to Noise 
		h 	  = gw.Signal(gw_prms, dets, approx)#*ExpT

		SNR2  = gw.ScalarProduct(det_b['freq'], det_b['Sn'],h,h)
		rho2 += SNR2
		GwData.append(h) ; SNR.append(np.sqrt(SNR2))
		#------------------------------------------------------
		# Computing Fisher Matrix
		#------------------------------------------------------

		if(not hide_info): print("\n\n\t Computing Fisher %s(%d) ..." % (det_b['name'],num_det))
		for i in range(Np):
			for j in range(i+1):
				if(not hide_info): print('Fisher: [%d-%d]' % (i+1,j+1), end='\r')
				xj = FreeParams[j]
				xi = FreeParams[i]
				Fij = gw.Fisher_ij(xi, xj, gw_prms, dets, approx, step_size, diff_order)
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
				value = gw.Doublet12(*[FreeParams[ii] for ii in Idx2],gw_prms,dets,approx, step_size, diff_order)
				for p in pmt([j,k]):
					jj, kk = p
					Db12[i][jj][kk] += value
			# --------------Doublet(2,2)--------------
			for Idx2 in indepD_22:
				i, j, k, l = Idx2
				value = gw.Doublet22(*[FreeParams[ii] for ii in Idx2],gw_prms,dets,approx, step_size, diff_order)
				for p1 in pmt([i,j]):
					ii, jj = p1
					for p2 in pmt([k,l]):
						kk, ll = p2
						Db22[ii][jj][kk][ll] += value
						if(not hide_info): print('(%s) [Sym] Doublet: '%det_b['name'], i, j, k, l,'...'*10, end='\r')
			#*****************************#*****************************#*****************************
			if(dali_method == 'Triplet'):
				indepT_13, indepT_23, indepT_33 = idxs_trip
				# --------------Triplet(1,3)--------------
				for Idx3 in indepT_13:
					i, j, k, l = Idx3
					value = gw.Triplet13(*[FreeParams[ii] for ii in Idx3],gw_prms,dets,approx, step_size, diff_order)
					for p in pmt([j,k,l]):
						jj, kk, ll = p
						Tp13[i][jj][kk][ll] += value
				# --------------Triplet(2,3)--------------
				for Idx3 in indepT_23:
					i, j, k, l, m = Idx3
					value = gw.Triplet23(*[FreeParams[ii] for ii in Idx3],gw_prms,dets,approx, step_size, diff_order)
					for p1 in pmt([i,j]):
						ii, jj = p1
						for p2 in pmt([k,l,m]):
							kk, ll, mm = p2
							Tp23[ii][jj][kk][ll][mm] += value
				# --------------Triplet(3,3)--------------
				for Idx3 in indepT_33:
					i, j, k, l, m, n = Idx3
					value = gw.Triplet33(*[FreeParams[ii] for ii in Idx3],gw_prms,dets,approx, step_size, diff_order)
					for p1 in pmt([i,j,k]):
						ii, jj, kk = p1
						for p2 in pmt([l,m,n]):
							ll, mm, nn = p2
							Tp33[ii][jj][kk][ll][mm][nn] += value
							if(not hide_info): print('(%s) [Sym] Triplet: '%det_b['name'], i, j, k, l, m, n,'...'*10, end='\r')
		# --------------# --------------# --------------					

	SNR_sum = np.sqrt(rho2)
	if(not hide_info):
		print(" >> SNR = ", SNR_sum)
		print(" >> 100/SNR = %.2e" % (100/SNR_sum) + '%')

	Doublet = [Db12, Db22]
	Triplet = [Tp13, Tp23, Tp33]

	return [SNR, GwData, Fisher , Doublet, Triplet]