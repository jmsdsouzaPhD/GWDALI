import sys, os
import GWDALI as gw
import numpy as np
import pandas as pd
import jax.numpy as jnp
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from time import time as now
from tqdm import trange

np.random.seed(0)

# Cosmic Explorer Coordinates
# https://cosmicexplorer.org/celocations.html
lon_CE = -125.0891
lat_CE = 45.64555

# Einstein Telecope Coordinates (Vaalserberg)
lon_ET = -(6. + 1./60 + 14./3600)
lat_ET = 50. + 45./60 + 16./3600

det_CE = {'name':'CE','lon':lon_CE, 'lat':lat_CE, 'rot':0, 'shape':90}
det_ET1 = {'name':'ET','lon':lon_ET, 'lat':lat_ET, 'rot':0, 'shape':60}
det_ET2 = {'name':'ET','lon':lon_ET, 'lat':lat_ET, 'rot':120, 'shape':60}
det_ET3 = {'name':'ET','lon':lon_ET, 'lat':lat_ET, 'rot':-120, 'shape':60}

detectors = [det_CE,det_ET1,det_ET2,det_ET3]

if False:
	xx, yy, F = gw.get_map(detectors,plot_map=False)

	rad = np.pi/180

	norm = mcolors.Normalize(vmin=0,vmax=np.max(F))

	fig = plt.figure()
	plt.subplot(111,projection='aitoff')
	plt.suptitle(f"min(F) = {np.min(F):.2f} ; max(F) = {np.max(F):.2f}")

	c = plt.contourf(xx*rad,yy*rad,F,cmap='magma',levels=np.linspace(0,np.max(F),21),norm=norm)
	plt.plot(lon_CE*rad, lat_CE*rad,'wo',mec='k')
	plt.plot(lon_ET*rad, lat_ET*rad,'wo',mec='k')

	cb = plt.colorbar(c,orientation='horizontal',label=' $F = \\sqrt{\\sum_a (F_{+,a}^2+F_{\\times,a}^2)}$')
	cb.set_ticks(np.linspace(0,np.max(F),6))

	plt.grid(ls='--',color='k',alpha=0.3)

	plt.tight_layout()
	fig.savefig('outputs/fig_patterns_ET+CE_flat.png')

#===================#========================

fig = gw.get_map(detectors,plot_map=True)
fig.savefig('outputs/fig_patterns_ET+CE_map.png')

plt.show()