import GWDALI as gw
import numpy as np

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

fig = gw.get_map(detectors,plot_map=True)
fig.savefig('outputs/fig_patterns_ET+CE_map.png')

plt.show()