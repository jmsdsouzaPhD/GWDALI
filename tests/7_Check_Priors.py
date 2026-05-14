import numpy as np
import GWDALI as gw
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(70,0.3)
dL = cosmo.luminosity_distance([.1,.5,1.]).value/1.e3

FreeParams = "dL,inv_dL,inv_dL2,inv_sqrtdL,inv_lnDL,lnDL,iota,cos_iota,psi,phi_coal,RA,Dec,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,chi_s,chi_a,m1,m2,M,q,inv_eta,ln_Mc,ln_eta,deltaM".split(',')

dict_dL = {"dL":dL, "inv_dL":1/dL, "inv_dL2":1/dL**2, "inv_sqrtdL":1/dL**.5}

Priors = gw.Priors(FreeParams,name=None,new_priors=None,plot=False)

fig = plt.figure(figsize=(12,8))
plt.suptitle("Default Priors",weight="bold")

for n, key in enumerate(FreeParams):
	plt.subplot(5,7,n+1)

	p = Priors[key]
	l1, l2 = p.minimum, p.maximum
	
	if key in dict_dL.keys():
		x = np.logspace(np.log10(l1),np.log10(l2),1000)
		y = p.prob(x)
		for u in dict_dL[key]:
			#plt.axvline(u,color='r',ls='--')
			plt.xscale('log')
	else:
		x = np.linspace(l1,l2,1000)
		y = p.prob(x)
	y/=max(y)
	
	#if key in "q,eta,inv_eta,deltaM".split(','):
	#	plt.xscale('log')
	
	
	plt.plot(x,y,'k-',lw=2)
	if min(y)!=max(y):
		l1 = max(y)/100
		xs = x[y>l1]
		plt.xlim(min(xs),max(xs))
		plt.xticks(rotation=45)

	plt.grid(alpha=.3)
	plt.xlabel(key)

plt.tight_layout()
fig.savefig("outputs/priors.png")
fig.savefig("outputs/priors.pdf")
plt.show()
print("OK!")