import numpy as np
import GWDALI as gw
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(70,0.3)
dL = cosmo.luminosity_distance([.1,.5,1.]).value/1.e3

FreeParams = "dL,inv_dL,inv_dL2,inv_sqrtdL,inv_lnDL,lnDL,iota,cos_iota,psi,phi_coal,RA,Dec,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,chi_s,chi_a,m1,m2,M,q,inv_eta,ln_Mc,ln_eta,deltaM".split(',')

keys="d_L,1/d_L,1/d_L^2,1/\\sqrt{d_L},1/ln(d_L),ln(d_L),\\iota,cos(\\iota),\\psi,\\phi_{coal},RA,Dec,t_{coal},M_c,\\eta,S_{x_1},S_{y_1},S_{z_1},S_{x_2},S_{y_2},S_{z_2},"
keys+="\\chi_s,\\chi_a,m_1,m_2,M,q,1/\\eta,ln(M_c),ln(\\eta),\\delta_M"
keys = [f"${x}$" for x in keys.split(',')]

TexKeys = {}
for i, key in enumerate(FreeParams):
	print(key,keys[i])
	TexKeys[key] = keys[i]

dict_dL = {"dL":dL, "inv_dL":1/dL, "inv_dL2":1/dL**2, "inv_sqrtdL":1/dL**.5}

Priors = gw.Priors(FreeParams,name=None,new_priors=None,plot=False)

fig = plt.figure(figsize=(16,10))
plt.suptitle("Default Priors",weight="bold")

for n, key in enumerate(FreeParams):
	plt.subplot(5,7,n+1)

	p = Priors[key]
	l1, l2 = p.minimum, p.maximum
	
	if key in dict_dL.keys():
		x = np.logspace(np.log10(l1),np.log10(l2),10000)
		y = p.prob(x)
	else:
		x = np.linspace(l1,l2,1000)
		y = p.prob(x)
	y/=max(y)
	
	if key in "dL,inv_dL,inv_dL2".split(','): plt.xscale('log')
	
	n = len(x)
	plt.plot(x,y,'k-',lw=2)
	#plt.fill_between(x,np,max(y)/100*np.ones(n),y,color='k',alpha=.3)
	if min(y)!=max(y):
		l1 = max(y)/100
		xs = x[y>l1]
		plt.xlim(min(xs),max(xs))
		plt.xticks(rotation=45)

	plt.grid(alpha=.3)
	plt.xlabel(TexKeys[key])

plt.tight_layout()
fig.savefig("outputs/priors.png")
fig.savefig("outputs/priors.pdf")
plt.show()
print("OK!")