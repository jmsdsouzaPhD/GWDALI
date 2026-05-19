import GWDALI as gw
import matplotlib.pyplot as plt

det1 = {'name':'det. 1','lon':-83, 'lat':30, 'rot':0, 'shape':90}
det2 = {'name':'det. 2','lon':-50, 'lat':-20, 'rot':30, 'shape':40}
det3 = {'name':'det. 3','lon':-65, 'lat':-35, 'rot':60, 'shape':70}
det4 = {'name':'det. 4','lon':10, 'lat':50, 'rot':90, 'shape':60}
det5 = {'name':'det. 5','lon':130, 'lat':-20, 'rot':120, 'shape':120}
det6 = {'name':'det. 6','lon':80, 'lat':25, 'rot':150, 'shape':50}

detectors = [det1,det2,det3,det4,det5,det6]

fig = gw.draw_detectors(detectors)
fig.savefig('outputs/fig_draw_detectors.png')

plt.show()