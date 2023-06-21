import numpy as np
from itertools import permutations

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