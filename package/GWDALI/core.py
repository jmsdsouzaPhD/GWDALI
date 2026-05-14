'''
BSD 3-Clause License

Copyright (c) 2026, Josiel Mendonça Soares de Souza

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import jit, vmap, config
from .lib.corner_plot import corner_plot
from .lib.Likelihoods import get_samples, Check_Priors
from .lib import AngTransf as geo
from .lib import Waveforms as wf
from .lib import Tensors as Tns
import matplotlib.pyplot as plt
from time import time as now
from tqdm import trange

import warnings, os, sys
import numpy as np
import h5py

warnings.filterwarnings("ignore")

PI = 3.141592653589793238462643383279502884
gE = 0.577215664901532860606512090082402431
Msun = 1.988409870698050731911960804878414216e30
pc = 3.085677581491367278913937957796471611e16
G = 6.67430e-11
c = 299792458.
Gpc = pc*1.e9

R_earth_equatiorial = 6.3781e6
R_earth_polar = 6.3568e6
R_earth = 6.378e6 # meters

A0 = 7.806521525937888e-23
# A0 := sqrt(5/24)*(G*Msun/c^3)^(5./6) / (pi^(2/3) * Gpc / c)

def GWDALI( GwPrms,
			detectors,
			FreeParams,
			new_priors = None,
			approx="TaylorF2",
			method="Doublet",
			sampler="nestle",
			diff_method="autodiff",
			dali_tensors=None, # It has to be an array [ Fisher, [Db12,Db22], [Tp13,Tp23,Tp33] ]
			step_size = [1.e-6],
			run_sampler = True,
			hide_info = False,
			npoints=300,
			nwalkers=None,
			ntemps=10,
			nburn=0.,
			limits=None,
			pos0=None,
			npool=None,
			verbose=True,
			remove_out=True,
			output_name=None,
			save_bilby_path=True,
			bilby_path="outputs_bilby/",
			enable_jax_waveforms=True,
			**kwargs):
		
	disable_jit = kwargs.get("disable_jit", False)
	jax.config.update("jax_disable_jit",disable_jit)
	if disable_jit: print(">> Disabling jax.jit()")

	#===================================================================================
	#============================== PRINT GWDALI_v1 SETTINGS ===============================
	#===================================================================================

	print("#"+"==="*10)
	print("... Running GWDALI_v1 ...")
	print("#"+"==="*10,"\n")
	print("\n Free Parameters: ", FreeParams)
	print(" Method: ", method)
	print(" Approx: ", approx,"\n")
	
	theta_keys, gw_values = zip(*GwPrms.items())

	ndim = len(FreeParams)
	print(">> Injections:")
	for key in theta_keys:
		print(f"\t {key}: {GwPrms[key]:.2f}")
	print("\n")
	#===================================================================================
	#============================ GET GW_DATA & COMPUTE SNR ============================
	#===================================================================================

	Data = []

	rho2 = 0. ; snr2 = 0.
	cont = 0

	det_conf_ref = [detectors[0][x] for x in "lon,lat,rot,shape".split(',')]
	for det in detectors:
		cont += 1
		Sn, freq = Tns.get_Sn(det['name'],**kwargs)
		det_conf = [ [det[x] for x in "lon,lat,rot,shape".split(',')] , det_conf_ref ]

		if enable_jax_waveforms:
			Gw_Signal = wf.load_waveforms(theta_keys,approx,**kwargs)[0]
			args = [*gw_values, det_conf, freq]
			h  = Gw_Signal[approx](gw_values, det_conf, freq)
		else:
			h_func = wf.build_waveform_strain_lal(theta_keys,approx,**kwargs)
			h = h_func(gw_values, det_conf,freq,approx,**kwargs)

		rho2 = Tns.ScalarProd(h,h,Sn,freq)
		snr2 += rho2
		print(det["name"],f"({cont}) SNR = {np.sqrt(rho2)}")
		Data.append(h)

	SNR = np.sqrt(snr2)

	#===================================================================================
	#============================== COMPUTE DALI TENSORS ===============================
	#===================================================================================

	Tensors, Results, Fisher = [], [], []
	Tensors_flat = []
	print("Computing Tensors ...")
	time_dali = None
	if(method in ["Fisher","Doublet","Triplet"]):
		if dali_tensors is None:
			print("Loading JAX derivatives ...")
			time_init = now()
			Tensors, time_dets_tns = get_dali_tensors(GwPrms,detectors,FreeParams,method,approx,enable_jax_waveforms,diff_method,step_size,hide_info,**kwargs)
			Fisher, Doublet, Triplet = Tensors.values()
			time_dali = now()-time_init
		else:
			print("\n<< DALI tensors loaded by the user >>\n")
			Fisher  = np.array( dali_tensors[0] )
			Doublet = [ np.array(x) for x in dali_tensors[1] ]
			Triplet = [ np.array(x) for x in dali_tensors[2] ]
			Tensors = [Fisher, Doublet, Triplet]
			print("tensors loaded successfully!")
		
		Fisher_flat  = np.ravel(Fisher)
		Doublet_flat = [ np.ravel(db) for db in Doublet ]
		Triplet_flat = [ np.ravel(tp) for tp in Triplet ]
		Tensors_flat = [Fisher_flat, Doublet_flat, Triplet_flat]
	print("\n>> DALI-Tensors OK!")

	#===================================================================================
	#=================================== RUN SAMPLER ===================================
	#===================================================================================

	Truths = [ GwPrms[fp] for fp in FreeParams]
	time_mcmc = None
	if(run_sampler):
		if nwalkers==None: nwalkers = 3*len(FreeParams)
		t_init = now()
		Results = get_samples(Data,Tensors_flat,detectors,GwPrms,FreeParams,\
							  approx,method,sampler,npoints,nwalkers,ntemps,nburn,verbose,\
							  new_priors,limits,pos0,npool,remove_out,output_name,save_bilby_path,
							  enable_jax_waveforms,bilby_path,**kwargs)
		time_mcmc = now() - t_init
	return Results, Truths, Tensors, Fisher, [time_dali, time_mcmc]

#=============================================#=============================================#=============================================#
#=============================================#=============================================#=============================================#

def get_gwjax_wf_functions():
	Gw_Signal, Gw_hphx, Gw_Signal_Real, Gw_Signal_Imag = wf.load_waveforms()
	return Gw_Signal, Gw_hphx, Gw_Signal_Real, Gw_Signal_Imag

def get_hphx(detectors,GwPrms,freq,approx, enable_jax_waveforms=True,**kwargs):
	if 'disable_jit' in kwargs.keys():
		if kwargs['disable_jit'] == True:
			jax.config.update("jax_disable_jit",True)
			print(">> Disabling jax.jit()")

	theta_keys, prms = zip( *GwPrms.items() )

	if enable_jax_waveforms:
		Gw_Signal, Gw_hphx, Gw_Signal_Real, Gw_Signal_Imag = wf.load_waveforms(theta_keys,approx)

	hphx_vec = []
	for det in detectors:
		args = list( GwPrms.values() )
		if enable_jax_waveforms:
			hp, hx = Gw_hphx[approx](args,freq)
		else:
			hphx_func = wf.build_waveform_hphx_lal(theta_keys,approx,**kwargs)
			hp, hx = hphx_func(prms,freq,approx,**kwargs)
		hphx_vec.append([hp,hx])
	return hphx_vec

def get_strain(detectors,GwPrms,freq,approx,enable_jax_waveforms=True,**kwargs):
	if 'disable_jit' in kwargs.keys():
		if kwargs['disable_jit'] == True:
			jax.config.update("jax_disable_jit",True)
			print(">> Disabling jax.jit()")

	theta_keys, gw_values = zip(*GwPrms.items())
	if enable_jax_waveforms:
		Gw_Signal, Gw_hphx, Gw_Signal_Real, Gw_Signal_Imag = wf.load_waveforms(theta_keys,approx,**kwargs)
	
	h_vec = []
	hr_vec = []
	hi_vec = []
	hide_cond = True
	if 'hide_info' in kwargs.keys():
		hide_cond = kwargs['hide_info']
	det_conf_ref = [detectors[0][x] for x in "lon,lat,rot,shape".split(',')]
	for det in detectors:
		if not hide_cond: print(f"[GWDALI_v1] get_strain >> det = {det['name']}")
		det_conf = [ [det[x] for x in "lon,lat,rot,shape".split(',')] , det_conf_ref ]

		if enable_jax_waveforms:
			args = [*gw_values, det_conf, freq]
			h  = Gw_Signal[approx](gw_values, det_conf, freq)
			hr = Gw_Signal_Real[approx](*args)
			hi = Gw_Signal_Imag[approx](*args)
		else:
			h_func = wf.build_waveform_strain_lal(theta_keys,approx,**kwargs)
			h = h_func(gw_values, det_conf,freq,approx,**kwargs)
			hr = h.real
			hi = h.imag
		if not hide_cond: print("\t\t >> h computed!")

		h_vec.append(h)
		hr_vec.append(hr)
		hi_vec.append(hi)
	return h_vec, hr_vec, hi_vec

def get_SNR(detectors,GwPrms,approx,enable_jax_waveforms=True,**kwargs):
	if 'disable_jit' in kwargs.keys():
		if kwargs['disable_jit'] == True:
			jax.config.update("jax_disable_jit",True)

	SNR_vec = [] ; snr2 = 0
	theta_keys, gw_values = zip(*GwPrms.items())
	if enable_jax_waveforms:
		Gw_Signal, Gw_hphx, Gw_Signal_Real, Gw_Signal_Imag = wf.load_waveforms(theta_keys,approx,**kwargs)

	det_conf_ref = [detectors[0][x] for x in "lon,lat,rot,shape".split(',')]
	for det in detectors:
		Sn, freq = Tns.get_Sn(det['name'],**kwargs)
		det_conf = [ [det[x] for x in "lon,lat,rot,shape".split(',')] , det_conf_ref ]
		
		if enable_jax_waveforms:
			args = [*gw_values, det_conf, freq]
			h  = Gw_Signal[approx](gw_values, det_conf, freq)
		else:
			h_func = wf.build_waveform_strain_lal(theta_keys,approx,**kwargs)
			h_lal = h_func(gw_values, det_conf,freq,approx,**kwargs)
			h = jnp.array( h_lal )

		rho2 = Tns.ScalarProd(h,h,Sn,freq)
		snr2 += rho2
		SNR_vec.append(np.sqrt(rho2))
	SNR_tot = np.sqrt(snr2)
	return SNR_vec, SNR_tot

def get_derivatives(FreeParams,approx,GwPrms,dets,freq,diff_order="first",diff_method="numdiff",
					step_size=1.e-6,full_tensor=True,enable_jax_waveforms=True,**kwargs):

	if 'disable_jit' in kwargs.keys():
		if kwargs['disable_jit'] == True:
			jax.config.update("jax_disable_jit",True)
			print(">> Disabling jax.jit()")
			
	jitgrad = kwargs.get("jitgrad", False)
	EarthRotation = kwargs.get("EarthRotation",False)

	gwkeys = list( GwPrms.keys() ) ; ndim = len(FreeParams)
	print(">> Computing Derivatives with respect to:",FreeParams,"\n")

	det_a = [dets[0][x] for x in "lon,lat,rot,shape".split(',')]
	if len(dets)>1:
		det_b = [dets[1][x] for x in "lon,lat,rot,shape".split(',')]
	else:
		det_b = det_a
	det_conf = [det_a, det_b]

	if diff_method=="autodiff" and enable_jax_waveforms: 
		print("\t loading <autodiff> derivatives ...")
		grads 	  = Tns.load_derivatives(gwkeys,approx,jitgrad=jitgrad,EarthRotation=EarthRotation)
		aux_diff  = [grads, det_conf, freq]
		print("\t <autodiff> derivatives loaded!\n")
	else: # numdiff
		aux_diff = [gwkeys,step_size, approx, det_conf, freq, enable_jax_waveforms]

	list_prms = list( GwPrms.values() )
	Diff_values, time_diff = Tns.Get_Derivatives(diff_order,diff_method,FreeParams,list_prms,aux_diff,freq,full_tensor=full_tensor,**kwargs)
	# time_vec = [ time_diff1[0],...,time_diff1[n], time_diff2[j], ..., time_diff3[N] ]
	# time_diff = [time_vec,dt_total]
	# Diff_values = Derivatives, [time_vec,dt_total]
	return Diff_values, time_diff

def get_dali_tensors(GwPrms,detectors,FreeParams,method,approx,enable_jax_waveforms=True,
				 diff_method="autodiff",step_size=[1.e-6,1.e-4,1.e-2],hide_info=False,**kwargs):
	
	step_size = list(step_size)
	if len(step_size)==0:
		print("\n>> Error: 'step_size' argument is empty.\n\tIt should be non-empty.")
		quit()
	elif len(step_size)==1:
		step_size = step_size*3
	elif len(step_size)==2:
		step_size = step_size+[step_size[-1]]
	else:
		step_size = step_size[:3]

	if 'disable_jit' in kwargs.keys():
		if kwargs['disable_jit'] == True:
			jax.config.update("jax_disable_jit",True)
			print(">> Disabling jax.jit() on gwjax.dali_tensors()")

	Fisher, Doublet, Triplet, time_dets_tns = Tns.get_tensors(GwPrms,approx,detectors,FreeParams,method,diff_method,step_size,enable_jax_waveforms,hide_info,**kwargs)
	Tensors = {"Fisher":Fisher, "Doublet":Doublet, "Triplet":Triplet}
	return Tensors, time_dets_tns

def Priors(FreeParams,name=None,new_priors=None,plot=False):
	return Check_Priors(FreeParams,name,new_priors,plot)

#=============================================#=============================================#=============================================#
#=============================================#=============================================#=============================================#

import matplotlib.image as mpimg
from pathlib import Path
path_map = Path(__file__).parent

def get_map(detectors,plot_map=False):
	print(">> Getting Detector Response ...")
	lons = np.linspace(-180,180,50)
	lats = np.linspace(-90,90,50)
	xx, yy = np.meshgrid(lons,lats)
	F2 = xx*0
	for det in detectors:
		for i in trange(50):
			for j in range(50):
				alpha_det, beta_det, psi_det = geo.AngTransf(0,0,lons[j],lats[i],det['lon'],det['lat'],det['rot'])
				Fp, Fx = wf.PatternFunc(alpha_det, beta_det, psi_det, det['shape'])
				Fp = np.array(Fp); Fx = np.array(Fx)
				F2[i][j] += Fp**2 + Fx**2
	F = np.sqrt(F2)

	if plot_map:
		img = mpimg.imread(path_map / "map.png")
		fig = draw_detectors(detectors)
		ax = plt.subplot(111)
		ax.set_xlim(-180, 180)
		ax.set_ylim(-90, 90)
		ax.set_xlabel("Lon")
		ax.set_ylabel("Lat")
		plt.grid(alpha=.3,ls='--')
		C = plt.contourf(xx,yy,F,levels=50,cmap='jet',zorder=1)
		ax.imshow(img,extent=[-180, 180, -90, 90],origin="upper",zorder=2)
		plt.colorbar(C,label="Detector Patterns Response")
		plt.tight_layout()
		return fig
	else:
		return xx, yy, F

def draw_detectors(dets):
	img = mpimg.imread(path_map / "map.png")
	fig = plt.figure(figsize=(12,6))
	ax = plt.subplot(111)
	ax.imshow(img,extent=[-180, 180, -90, 90],origin="upper")

	for det in dets:
		x0, y0 = det['lon'],det['lat']

		r = (det['rot']-90)*np.pi/180
		sh = det['shape']*np.pi/180
		
		ux = (np.pi/4 - sh/2)
		uy = (np.pi/4 + sh/2)

		l = 10.
		x = [x0+l*np.cos(r + ux),x0,x0+l*np.cos(r + uy)]
		y = [y0+l*np.sin(r + ux),y0,y0+l*np.sin(r + uy)]

		ax.plot(x,y,'o-',mec='k',lw=3,label=det['name'],zorder=3)

	plt.legend()
	ax.set_xlim(-180, 180)
	ax.set_ylim(-90, 90)
	ax.set_xlabel("Lon")
	ax.set_ylabel("Lat")
	plt.grid(alpha=.3,ls='--')
	plt.tight_layout()
	return fig