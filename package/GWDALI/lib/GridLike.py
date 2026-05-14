import numpy as np
from time import time as now

def Call_Grid(likelihood,priors,npoints,limits=None,scale='linear'):
	keys = list(priors.keys())
	ndim = len(keys)
	Prms = {}

	#if limits is None:
	#	limits = {}
	#	for key in keys:
	#		limits[key] = priors[key].minimum , priors[key].maximum
	if limits is None:
		limits = {}
		for key in keys:
			limits[key] = None
	else:
		for key in keys:
			if not key in limits.keys():
				limits[key] = None

	shape = (npoints,)*ndim
	Pr, Theta = [], []
	for key in keys:
		if limits[key] is None:
			x1, x2 = priors[key].minimum , priors[key].maximum
			Prms[key] = np.linspace(x1,x2,npoints)
		else:
			Prms[key] = limits[key]
		Pr.append( priors[key].prob(Prms[key]) )
		Theta.append(Prms[key])

	Prior = Pr[0]
	for i in range(1,len(keys)):
		Prior = np.tensordot(Prior, Pr[i], axes=0)
	log_prior = np.log(Prior)

	logP = np.zeros(shape)
	all_idxs = list(np.ndindex(shape))
	N = len(all_idxs)
	t_init = now() ; likelihood_times = []
	for cont, idxs in enumerate(all_idxs):
		prms = [ Theta[i][j] for i, j in enumerate(idxs) ]
		likelihood.parameters = {}
		for n, key in enumerate(keys):
			likelihood.parameters[key] = prms[n]
		t0 = now()
		logP[tuple(idxs)] = likelihood.log_likelihood() + log_prior[tuple(idxs)]
		dT = now() - t_init ; likelihood_times.append(now()-t0)
		dt = dT/(cont+1)

		T = dT 
		t = (N-cont)*dt

		print(f"Grid: ({cont+1}:{N} , {int(100*cont/N)}%)  [{int(T//60)}min:{int(T%60)}sec - {int(t//60)}min:{int(t%60)}sec]", end='\r')
	logN = np.max(logP)
	Posterior = np.exp( logP.T - logN )
	return Prms, Posterior, np.array(likelihood_times)