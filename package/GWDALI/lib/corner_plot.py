import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter as gf

def corner_plot(Samples,labels=None,truths=None,CLs=None,bins=30,truth_color='r',smooth=1,smooth1d=1,cmap='jet'):
	ndim = len(Samples.T)
	fig = plt.figure(figsize=(6,6))
	#print("figsize = ", ndim*2,ndim*2)
	for i in range(ndim):
		xi = Samples[:,i]
		ax = plt.subplot(ndim,ndim,i*ndim+i+1)
		#plt.hist(xi,bins=bins,histtype='step',lw=2,color='k')
		u, v = np.histogram(xi,bins=bins)
		dv = v[1]-v[0] ; v = v[1:]-dv/2
		u = u/max(u)
		y = gf(u,smooth1d) ; y/=np.max(y)
		plt.plot(v,y,'k-',lw=2)
		plt.fill_between(v,np.zeros(len(v)),y,color='k',alpha=0.3)
		if(truths!=None): ax.vlines(truths[i],0,1.1,color=truth_color)
		ax.set_ylim(0,1.1)
		ax.set_xlim(min(xi),max(xi))
		ax.grid(ls='--',alpha=0.3)
		ax.set_yticklabels([])
		if(all([s != None for s in [labels,CLs]])):
			CLs = sorted(CLs)
			err_p = [0.5+cl/2 for cl in CLs]
			err_m = [0.5-cl/2 for cl in CLs]
			qs = sorted(err_m + [0.5] + err_p)
			qs = np.quantile(xi,qs) ; n = int(len(qs)/2)
			err_u = qs[n+1]-qs[n]
			err_d = qs[n]-qs[n-1]
			ax.vlines(qs[n-1],0,1.1,color='k',ls='--',lw=1)
			ax.vlines(qs[n],0,1.1,color='k',ls='--',lw=1)
			ax.vlines(qs[n+1],0,1.1,color='k',ls='--',lw=1)
			plt.title('%s=$%.2f^{+%.2f}_{-%.2f}$' % (labels[i],qs[n],err_u,err_d))
			if(i==ndim-1): ax.set_xlabel(labels[i])
		for j in range(i):
			xj = Samples[:,j]
			ax = plt.subplot(ndim,ndim,i*ndim+j+1)
			H, xe, ye = np.histogram2d(xj,xi,bins=bins)
			dxe = xe[1]-xe[0] ; dye = ye[1]-ye[0]
			xe = xe[1:]-dxe/2 ; ye = ye[1:]-dye/2
			xx, yy = np.meshgrid(xe,ye)

			plt.imshow(H.T,origin='lower',cmap=cmap,extent=(min(xj),max(xj),min(xi),max(xi)),aspect='auto',interpolation='bicubic')
			plt.plot(xj,xi,color='w',marker='.',lw=0,ms=.1,alpha=.5,zorder=1)

			H_vec = np.sort(np.ravel(H))
			norm = np.sum(H_vec)
			H_vec = H_vec/norm ; H = H/norm
			Cum = np.array([np.sum(H_vec[k:]) for k in range(len(H_vec))])
			func1 = interp1d(Cum,H_vec,bounds_error=False,fill_value='extrapolate')

			vals = func1(np.flip(CLs))
			try: plt.contour(xx,yy,gf(H.T,smooth),levels=vals,colors='w',zorder=2,linewidths=1)
			except: pass
			plt.grid(ls='--',alpha=.3)
			ax.set_xlim(min(xj),max(xj))
			ax.set_ylim(min(xi),max(xi))
			if(i!=(ndim-1)): ax.set_xticklabels([])
			else: ax.set_xlabel(labels[j])
			if(j!=0): ax.set_yticklabels([])
			else: ax.set_ylabel(labels[i])
			if(truths!=None):
				ax.vlines(truths[j],min(xi),max(xi),color=truth_color)
				ax.hlines(truths[i],min(xj),max(xj),color=truth_color)
				ax.plot(truths[j],truths[i],'s',color=truth_color,mec='k')
	fig.subplots_adjust(wspace=0.03,hspace=0.03)
	return fig