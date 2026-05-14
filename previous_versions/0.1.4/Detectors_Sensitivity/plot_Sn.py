import numpy as np
import matplotlib.pyplot as plt
import os

files = ['Sn_L.txt','Sn_V.txt','Sn_K.txt','Sn_ET.txt','Sn_CE.txt']
fig = plt.figure()

labels = {}
labels['Sn_L'] = 'aLIGO'
labels['Sn_V'] = 'aVirgo'
labels['Sn_K'] = 'KAGRA'
labels['Sn_ET'] = 'ET'
labels['Sn_CE'] = 'CE'

ls = {}
ls['Sn_L'] = '-'
ls['Sn_V'] = '--'
ls['Sn_K'] = ':'
ls['Sn_ET'] = '-'
ls['Sn_CE'] = '--'

colors = {}
colors['Sn_L'] = 'b'
colors['Sn_V'] = 'b'
colors['Sn_K'] = 'b'
colors['Sn_ET'] = 'r'
colors['Sn_CE'] = 'r'

for f in files:
	if(f.endswith('.txt')):
		data = np.loadtxt(f)
		x = data[:,0]
		y = np.sqrt(x)*data[:,1]
		ref = f.replace('.txt','')

		plt.loglog(x,y,color=colors[ref],ls=ls[ref],lw=3,label=labels[ref])
		print(f, len(x))

plt.legend()
plt.grid(which='both',alpha=0.3)
plt.xlabel('frequency [Hz]')
plt.ylabel('Characteristic Strain $\sqrt{f\cdot S_n(f)}$')

fig.savefig('Sn.png')
fig.savefig('Sn.pdf')

plt.show()

