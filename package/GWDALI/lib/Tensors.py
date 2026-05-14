import numpy as np
import jax
import jax.numpy as jnp
from jax import jit, grad, jacrev, vmap, config
from itertools import permutations
from time import time as now
from pathlib import Path
from tqdm import trange
import warnings, sys, os

#warnings.filterwarnings("ignore")

from . import Waveforms as wf
from . import Auxiliar as sym

Kappa = 7.806521525937888e-23
#Kappa := sqrt(5/24)*(G*Msun/c^3)^(5./6) / (pi^(2/3) * Gpc / c)

dets = ['ET','CE','aLIGO','aVirgo','KAGRA']
path = Path(__file__).parent / '../Sensitivities/'

Sensitivity = {}
Sensitivity['LIGO']  = np.loadtxt(path / 'Sn_L.txt')
Sensitivity['Virgo'] = np.loadtxt(path / 'Sn_V.txt')
Sensitivity['Kagra'] = np.loadtxt(path / 'Sn_K.txt')
Sensitivity['ET']    = np.loadtxt(path / 'Sn_ET.txt')
Sensitivity['CE']    = np.loadtxt(path / 'Sn_CE.txt')
Sensitivity['CE_40km']    = np.loadtxt(path / 'Sn_CE_40km.txt')
Sensitivity['CE_20km']    = np.loadtxt(path / 'Sn_CE_20km.txt')
Keys= Sensitivity.keys()

try:
	from jax.scipy.integrate import trapezoid
except:
	print("Fail on loading jax.scypi.integrate.trapezoid")
	@jit
	def trapezoid(y,x):
		dx = jnp.diff(x)
		return jnp.sum( 0.5 * (y[1:]+y[:-1])*dx )

def get_Sn(name,**kwargs):
	if name in Sensitivity.keys():
		freq = jnp.array( Sensitivity[name][:,0] )
		Sn   = jnp.array( Sensitivity[name][:,1]**2)
	else:
		Sn, freq = kwargs["new_detector"]
		Sensitivity[name] = np.transpose([jnp.array(freq),jnp.sqrt(Sn)])
	return Sn, freq

#--------------------------------------------------------------------
# -------------------------- GW Parameters --------------------------
#--------------------------------------------------------------------

Keys = "dL,iota,psi,phi_coal,RA,Dec,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,\
ln_dL,inv_dL,inv_dL2,inv_sqrtdL,inv_lndL,cos_iota,inv_eta,ln_eta,ln_Mc,\
m1,m2,M,q,deltaM,chi_s,chi_a".split(',')
'''
	-------std_args-------
	[0]  dL
	[1]  iota
	[2]  psi
	[3]  phi_coal
	[4]  RA
	[5]  Dec
	[6]  t_coal
	[7]  Mc
	[8]  eta
	[9]  sx1
	[10] sy1
	[11] sz1
	[12] sx2
	[13] sy2
	[14] sz2
	-------var_args-------
	[15] ln_dL
	[16] inv_dL
	[17] inv_dL2
	[18] inv_sqrtdL
	[19] inv_lndL
	[20] cos_iota
	[21] inv_eta
	[22] ln_eta
	[23] ln_Mc
	[24] m1
	[25] m2
	[26] M
	[27] q
	[28] deltaM : Mass Asymmetry Parameter
	[29] chi_s
	[30] chi_a
'''

#----------------------------
# 11 independent parameters
#----------------------------
# Number of derivatives:
# 	>>  11 1st derivatives
# 	>>  66 2nd derivatives
# 	>> 286 3rd derivatives

#=============================================================================
#-----------------------------Autodiff jax.grad()-----------------------------
#=============================================================================

def gradient(function,argnums,jitgrad):
	if jitgrad:
		return jax.jit( grad(function, argnums) )
	else:
		return grad(function, argnums)

def load_derivatives(theta_keys,approx,jitgrad,**kwargs):
	# Precompute first-order gradients
	print(f"\t\t>> jitgrad = {jitgrad}...")
	
	Gw_Signal_Real, Gw_Signal_Imag = wf.load_waveforms(theta_keys,approx,**kwargs)[2:]
	globals()["Gw_Signal_R"] = Gw_Signal_Real
	globals()["Gw_Signal_I"] = Gw_Signal_Imag
	#==============================
	# Using jax.grad()
	#==============================
	for ri in ["R","I"]:
		
		globals()[f"grad1_{ri}"] = {    key: gradient(globals()[f"Gw_Signal_{ri}"][approx], argnums=i, jitgrad=jitgrad )  for i, key in enumerate(theta_keys) }
		
		globals()[f"grad2_{ri}"] = {	key1: {key2: gradient(globals()[f"grad1_{ri}"][key1], argnums=j, jitgrad=jitgrad ) for j, key2 in enumerate(theta_keys)}
		    		for i, key1 in enumerate(theta_keys) }

		globals()[f"grad3_{ri}"] = {	key1: { key2: {key3: gradient(globals()[f"grad2_{ri}"][key1][key2], argnums=k, jitgrad=jitgrad ) for k, key3 in enumerate(theta_keys)}
			        		for j, key2 in enumerate(theta_keys) }
		    		for i, key1 in enumerate(theta_keys) }

	return [grad1_R, grad1_I], [grad2_R, grad2_I], [grad3_R, grad3_I]

# Function to compute first-order derivatives
def aDiff1(prms0, x, grads, det_conf, freq, **kwargs):
    grad1_R, grad1_I = grads[0]
    #params = [p[key] for key in Keys]
    diffR = vmap(lambda f: grad1_R[x](*prms0, det_conf, f))(freq)
    diffI = vmap(lambda f: grad1_I[x](*prms0, det_conf, f))(freq)
    return diffR + 1.j * diffI

# Function to compute second-order derivatives
def aDiff2(prms0, xi, xj, grads, det_conf, freq, **kwargs):
    grad2_R, grad2_I = grads[1]
    #params = [p[key] for key in Keys]
    diffR = vmap(lambda f: grad2_R[xi][xj](*prms0, det_conf, f))(freq)
    diffI = vmap(lambda f: grad2_I[xi][xj](*prms0, det_conf, f))(freq)
    return diffR + 1.j * diffI

# Function to compute third-order derivatives
def aDiff3(prms0, xi, xj, xk, grads, det_conf, freq, **kwargs):
    grad3_R, grad3_I = grads[2]
    #params = [p[key] for key in Keys]
    diffR = vmap(lambda f: grad3_R[xi][xj][xk](*prms0, det_conf, f))(freq)
    diffI = vmap(lambda f: grad3_I[xi][xj][xk](*prms0, det_conf, f))(freq)
    return diffR + 1.j * diffI

#=============================================================================
#-----------------------Derivative: Finite Differences-----------------------
#=============================================================================

lims = {}
lims["iota"] = 0, np.pi
lims["RA"] = -180, 180
lims["Dec"] = -90, 90
lims["eta"] = 0, 0.25
lims["sz1"] = -1, 1
lims["sz2"] = -1, 1

def split_prms(prms0,theta_keys,eps,x_key):
	idx = theta_keys.index(x_key)
	x0 = prms0[idx]
	dx = max([eps,eps*jnp.abs(x0)])
	prms1 = prms0.copy() ; prms2 = prms0.copy()
	
	x1 = x0-dx/2 ; x2 = x0 + dx/2
	if x_key in lims.keys():
		if x1<=lims[x_key][0]: # avoid x extrapolate parameters_limits
			x1 = lims[x_key][0] ; x2 = x0+dx
		elif(x2>=lims[x_key][1]): # avoid x extrapolate parameters_limits
			x2 = lims[x_key][1] ; x1 = x0-dx

	prms1[idx] = x1
	prms2[idx] = x2

	return prms1, prms2, dx

def numDiff1(prms0, x, *aux_diff_prms,**kwargs):
	theta_keys, eps, approx, det_conf, freq, enable_jax_waveforms = aux_diff_prms
	prms1, prms2, dx = split_prms(prms0,theta_keys,eps,x)
	if enable_jax_waveforms:
		#======================
		# >> jax waveforms <<
		#======================
		Gw_Signal = wf.load_waveforms(theta_keys,approx,**kwargs)[0]
		
		y2 = Gw_Signal[approx]( prms2,det_conf,freq )
		y1 = Gw_Signal[approx]( prms1,det_conf,freq )
	else: 
		#======================
		# >> lal waveforms <<
		#======================
		#args = [*[prms[key] for key in Keys],det_conf,freq]
		#h = wf.Signal_lal(*list(args+[approx]))
		args1 = [*prms1,det_conf,freq,approx]
		args2 = [*prms2,det_conf,freq,approx]
		#y2 = wf.Signal_lal(*args1,**kwargs)
		#y1 = wf.Signal_lal(*args2,**kwargs)
		h_func = wf.build_waveform_strain_lal(theta_keys,approx,**kwargs)
		y1 = h_func(prms1, det_conf,freq,approx,**kwargs)
		y2 = h_func(prms2, det_conf,freq,approx,**kwargs)

	return (y2-y1)/dx

def numDiff2(prms0, xi, xj, *aux_diff_prms,**kwargs):	
	theta_keys, eps, approx, det_conf, freq, enable_jax_waveforms = aux_diff_prms
	prms1, prms2, dx = split_prms(prms0,theta_keys,eps,xi)
	y2 = numDiff1(prms2, xj, *aux_diff_prms,**kwargs)
	y1 = numDiff1(prms1, xj, *aux_diff_prms,**kwargs)
	return (y2-y1)/dx

def numDiff3(prms0, xi, xj, xk, *aux_diff_prms,**kwargs):
	theta_keys, eps, approx, det_conf, freq, enable_jax_waveforms = aux_diff_prms
	prms1, prms2, dx = split_prms(prms0,theta_keys,eps,xi)
	y2 = numDiff2(prms2, xj, xk, *aux_diff_prms,**kwargs)
	y1 = numDiff2(prms1, xj, xk, *aux_diff_prms,**kwargs)
	return (y2-y1)/dx

#=============================================================================

Diff1 = {"numdiff":numDiff1 , "autodiff": aDiff1}
Diff2 = {"numdiff":numDiff2 , "autodiff": aDiff2}
Diff3 = {"numdiff":numDiff3 , "autodiff": aDiff3}

@jit
def ScalarProd(A,B,Sn,freq):
	return 4*jnp.real( trapezoid( A*jnp.conj(B)/Sn, freq ) )

#=============================================================================
#------------------------------Store Derivatives------------------------------
#=============================================================================

def Get_Derivatives(diff_order,diff_method,FreeParams,list_prms,aux_diff,freq,full_tensor=True,**kwargs):
	ndim = len(FreeParams)
	# Autodiff >> aux_diff = [grads, det_conf, freq] # where >> grads = load_derivatives(approx,jitgrad)
	# Numdiff  >> aux_diff = [step_size, approx, det_conf, freq, enable_jax_waveforms]

	time_vec = []
	if diff_order == "first":
		print(f"\n\n # Computing First-Order Derivatives (diff_method = {diff_method}) .....\n")
		t_init = now()
		#Diff1_values = {Xi: Diff1[diff_method](list_prms, Xi, *aux_diff) for Xi in FreeParams}
		Diff1_values = {}
		for n, Xi in enumerate(FreeParams):
			ti = now()
			Diff1_values[Xi] = Diff1[diff_method](list_prms, Xi, *aux_diff,**kwargs) #; quit()
			dt1 = now() - ti
			
			print(f"\t ({Xi})\t [{n+1}-{ndim}]\t Diff[{Xi}] OK! \t({int(now()-ti)} seconds)")
			time_vec.append(dt1)

		D1_values = list( Diff1_values.values() )
		Derivatives = D1_values.copy()
		dt_total = now() - t_init
		print(f"\n>> First Derivatives Concluded in {int(dt_total//60)}min:{int(dt_total%60)}sec!\n")

	elif diff_order == "second":
		print("Computing Second-Order Derivatives ...")
		t_init = now() ; cont2 = 0 ; Ntot2 = int(ndim*(ndim+1)/2)
		D2_values = np.zeros([ndim,ndim,len(freq)],dtype=complex)
		D2_vector = []
		for i in range(ndim):
			for j in range(i+1):
				Xi = FreeParams[i]
				Xj = FreeParams[j]

				ti = now() ; cont2 += 1
				val2 = Diff2[diff_method](list_prms,Xi,Xj,*aux_diff,**kwargs)
				D2_vector.append(val2)
				print(f"\t ({Xi},{Xj})\t [{cont2}-{Ntot2}]\t Diff[{Xi},{Xj}] OK! \t({int(now()-ti)} seconds)")
				dt2 = now() - ti
				time_vec.append(dt2)
				for ii, jj in sym.pmt([i,j]):
					D2_values[ii][jj] = val2
		if full_tensor:
			Derivatives = D2_values.copy()
		else:
			Derivatives = np.array(D2_vector)
		dt_total = now() - t_init
		print(f"\n>> Second Derivatives Concluded in {int(dt_total//60)}min:{int(dt_total%60)}sec!\n")

	elif diff_order == "third":		
		print("Computing Third-Order Derivatives...")
		t_init = now() ; cont3 = 0 ; Ntot3 = int(ndim*(ndim+1)*(ndim+2)/(2*3))
		D3_values = np.zeros([ndim,ndim,ndim,len(freq)],dtype=complex)
		D3_vector = []
		for i in range(ndim):
			for j in range(i+1):
				for k in range(j+1):
					Xi = FreeParams[i]
					Xj = FreeParams[j]
					Xk = FreeParams[k]
					
					ti = now() ; cont3 += 1
					val3 = Diff3[diff_method](list_prms,Xi,Xj,Xk,*aux_diff,**kwargs)
					D3_vector.append(val3)
					print(f"\t ({Xi},{Xj},{Xk})\t [{cont3}-{Ntot3}]\t Diff[{Xi},{Xj},{Xk}] OK! \t({int(now()-ti)} seconds)")
					dt3 = now()-ti
					time_vec.append(dt3)

					for ii, jj, kk in sym.pmt([i,j,k]):
						D3_values[ii][jj][kk] = val3
		if full_tensor:
			Derivatives = D3_values.copy()
		else:
			Derivatives = np.array(D3_vector)
		dt_total = now() - t_init
		print(f"\n>> Third Order Derivatives Concluded in {int(dt_total//60)}min:{int(dt_total%60)}sec!\n")
	else:
		print(f"\n Error: Invalid <diff_method>: {diff_method}\n")
		quit()

	return Derivatives, [time_vec,dt_total]

def get_tensors(gwprms,approx,detectors,FreeParams,method,diff_method,step_size,enable_jax_waveforms=True,hide_info=False,**kwargs):
	step_size1, step_size2, step_size3 = step_size

	jitgrad = kwargs.get("jitgrad", False)
	EarthRotation = kwargs.get("EarthRotation",False)

	if hide_info: sys.stdout = open(os.devnull, 'w')

	if 'save_tensors' in kwargs.keys():
		path, index = kwargs['save_tensors']
		os.makedirs(f'{path}',exist_ok=True)

	ndim = len(FreeParams)
	Fisher = np.zeros([ndim,ndim])
	Db12 = np.zeros([ndim,ndim,ndim])
	Db22 = np.zeros([ndim,ndim,ndim,ndim])
	Tp13 = np.zeros([ndim,ndim,ndim,ndim])
	Tp23 = np.zeros([ndim,ndim,ndim,ndim,ndim])
	Tp33 = np.zeros([ndim,ndim,ndim,ndim,ndim,ndim])
	
	Doublet = [Db12, Db22]
	Triplet = [Tp13, Tp23, Tp33]

	gwkeys, gw_prms = zip(*gwprms.items())
	list_prms = np.array(gw_prms)

	cont = 0

	time_det_tns = []
	det_b = [detectors[0][x] for x in "lon,lat,rot,shape".split(',')] 
	for det in detectors:
		dt_fisher = 0
		dt_db12 = 0
		dt_db22 = 0
		dt_tp13 = 0
		dt_tp23 = 0
		dt_tp33 = 0
		cont+=1
		try:
			Sn, freq = get_Sn(det['name'])
			print(f"\n*********************************** [{cont}] Detector: {det['name']} ***********************************\n")
		except:
			print("Invalid Detector Name!!!")
			quit()

		det_a = [det[x] for x in "lon,lat,rot,shape".split(',')]
		det_conf = [det_a, det_b] 

		if not enable_jax_waveforms:
			Sn = np.array(Sn)
			freq = np.array(freq)

		if diff_method=="autodiff" and enable_jax_waveforms: 
			print("\t loading <autodiff> derivatives ...")
			grads 	  = load_derivatives(gwkeys,approx,jitgrad=jitgrad,EarthRotation=EarthRotation)
			for i in range(3):
				globals()[f'aux_diff{i+1}']  = [grads, det_conf, freq]

			print("\t <autodiff> derivatives loaded!\n")
		else: # numdiff
			for i in range(3):
				dx = locals()[f'step_size{i+1}']
				globals()[f'aux_diff{i+1}'] = [gwkeys, dx, approx, det_conf, freq, enable_jax_waveforms]

		if(method in ["Fisher","Doublet","Triplet"]):								
			time_diff1 = [[0.], 0.]
			time_diff2 = [[0.], 0.]
			time_diff3 = [[0.], 0.]
			print(f"\n==================================== Computing DALI Tensors ({method}) ====================================\n")
			
			D1_values, time_diff1 = Get_Derivatives("first",diff_method,FreeParams,list_prms,aux_diff1,freq,**kwargs)
			time_f = now() ; print('\n')
			for i in range(ndim):
				for j in range(i+1):
					for ii, jj in sym.pmt([i,j]):
						D_i = D1_values[ii]
						D_j = D1_values[jj]
						Fisher[ii][jj] += ScalarProd(D_i,D_j,Sn,freq)

					print(f"\tFisher(1|1) [{i},{j}]","..."*10,end="\r")
			print('\n')
			dt_fisher = now()-time_f

			if 'save_tensors' in kwargs.keys():
				path, index = kwargs['save_tensors']
				np.savez_compressed(path+f'Fisher_det{cont}_{index}.npz',Fisher=Fisher,dt_fisher=np.array([time_diff1[1],dt_fisher]))


			if(method in ["Doublet","Triplet"]):
				D2_values, time_diff2 = Get_Derivatives("second",diff_method,FreeParams,list_prms,aux_diff2,freq,**kwargs)
				# Computing Doublet/Triplet Tensors ...
				
				idxs_doub, idxs_trip = sym.get_indep_indexies(ndim)
				indepD_12, indepD_22 = idxs_doub # Doublet indexies symetries

				# +++++++++++++++++++++++++++++++++++Doublet(2,1):+++++++++++++++++++++++++++++++++++
				time_db12 = now() ; print('\n')
				for Idx2 in indepD_12:
					i, j, k = Idx2
					D_i  = D1_values[i]
					D_jk = D2_values[j][k]
					val = ScalarProd(D_i,D_jk,Sn,freq)
					for jj, kk in sym.pmt([j,k]):
						Db12[i][jj][kk] += val
						print(f'Doublet(1|2) [{i,j,k}] ','...'*10, end='\r')
				dt_db12 = now()-time_db12
				# +++++++++++++++++++++++++++++++++++Doublet(2,2):+++++++++++++++++++++++++++++++++++
				time_db22 = now() ; print()
				for Idx2 in indepD_22:
					i, j, k, l = Idx2
					D_ij = D2_values[i][j]
					D_kl = D2_values[k][l]
					val = ScalarProd(D_ij,D_kl,Sn,freq)
					for ii, jj in sym.pmt([i,j]):
						for kk, ll in sym.pmt([k,l]):
							Db22[ii][jj][kk][ll] += val
							print(f'Doublet(2|2) [{i,j,k,l}] ','...'*10, end='\r')
				dt_db22 = now()-time_db22 ; print('\n')
				if 'save_tensors' in kwargs.keys():
					path, index = kwargs['save_tensors']
					np.savez_compressed(path+f'Doublet_det{cont}_{index}.npz',Db12=Db12,Db22=Db22,dt_doublet=np.array([time_diff2[1],dt_db12,dt_db22]))

			if method == "Triplet":

				D3_values, time_diff3 = Get_Derivatives("third",diff_method,FreeParams,list_prms,aux_diff3,freq,**kwargs)
				indepT_13, indepT_23, indepT_33 = idxs_trip
				# +++++++++++++++++++++++++++++++++++Triplet(1,3):+++++++++++++++++++++++++++++++++++
				time_tp13 = now()
				for Idx3 in indepT_13:
					i, j, k, l = Idx3
					D_i   = D1_values[i]
					D_jkl = D3_values[j][k][l]
					val = ScalarProd(D_i,D_jkl,Sn,freq)
					for jj, kk, ll in sym.pmt([j,k,l]):
						Tp13[i][jj][kk][ll] += val
						print(f'Triplet(1|3) [{i,j,k,l}] ','...'*10, end='\r')
				dt_tp13 = now() - time_tp13 ; print()
				# +++++++++++++++++++++++++++++++++++Triplet(2,3):+++++++++++++++++++++++++++++++++++
				time_tp23 = now()
				for Idx3 in indepT_23:
					i, j, k, l, m = Idx3
					D_ij  = D2_values[i][j]
					D_klm = D3_values[k][l][m]
					val = ScalarProd(D_ij,D_klm,Sn,freq)
					for ii, jj in sym.pmt([i,j]):
						for kk, ll, mm in sym.pmt([k,l,m]):
							Tp23[ii][jj][kk][ll][mm] += val
							print(f'Triplet(2|3) [{i,j,k,l,m}] ','...'*10, end='\r')
				dt_tp23 = now() - time_tp23 ; print()
				# +++++++++++++++++++++++++++++++++++Triplet(3,3):+++++++++++++++++++++++++++++++++++
				time_tp33 = now()
				for Idx3 in indepT_33:
					i, j, k, l, m,n = Idx3
					D_ijk = D3_values[i][j][k]
					D_lmn = D3_values[l][m][n]
					val = ScalarProd(D_ijk,D_lmn,Sn,freq)
					for ii, jj, kk in sym.pmt([i,j,k]):
						for ll, mm, nn in sym.pmt([l,m,n]):
							Tp33[ii][jj][kk][ll][mm][nn] += val
							print(f'Triplet(3|3) [{i,j,k,l,m,n}] ','...'*10, end='\r')
				dt_tp33 = now() - time_tp33 ; print()

				if 'save_tensors' in kwargs.keys():
					path, index = kwargs['save_tensors']
					np.savez_compressed(path+f'Triplet_det{cont}_{index}.npz',Tp13=Tp13,Tp23=Tp23,Tp33=Tp33,dt_triplet=np.array([time_diff3[1],dt_tp13,dt_tp23,dt_tp33]))
			
			time_det_tns.append([	
									[time_diff1,time_diff2,time_diff3],
									[dt_fisher,dt_db12,dt_db22,dt_tp13,dt_tp23,dt_tp33]
								])	
		else:
			print("Invalid Method!")
			quit()

	Doublet = [Db12, Db22]
	Triplet = [Tp13, Tp23, Tp33]
	
	return Fisher, Doublet, Triplet, time_det_tns