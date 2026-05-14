import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter as gf
from scipy.interpolate import interp1d

def Add_Ellipse(cov2, center, ax, fill=False):
	det = np.linalg.det(cov2)

	W = cov2[0,0]+cov2[1,1]
	a2 = 0.5*( W + np.sqrt(W*W - 4*det) )
	b2 = 0.5*( W - np.sqrt(W*W - 4*det) )

	a = np.sqrt(a2) ; b = np.sqrt(b2)
	sin_theta = ( cov2[0,0] - a2 )
	cos_theta = cov2[0,1]
	angle = -np.arctan2( sin_theta , cos_theta )

	def get_xy(t,a,b):
		xe = center[0] + a*np.cos(t)*np.cos(angle) - b*np.sin(t)*np.sin(angle)
		ye = center[1] + a*np.cos(t)*np.sin(angle) + b*np.sin(t)*np.cos(angle)	
		return xe, ye

	if(fill):
		t = np.linspace(0,2*np.pi,100)
		x1, y1 = get_xy(t,a,b)
		x2, y2 = get_xy(t,2*a,2*b)
		ax.plot(x1,y1,'k-',lw=2,zorder=3)
		ax.plot(x2,y2,'k-',lw=2,zorder=3)
	else:
		t = np.linspace(0,2*np.pi,100)
		x1, y1 = get_xy(t,a,b)
		x2, y2 = get_xy(t,2*a,2*b)
		ax.plot(x1,y1,'r--',zorder=3)
		ax.plot(x2,y2,'r--',zorder=3)

	pass

#--------------------------------------------
# Corner Plot
#--------------------------------------------

def Curves_CL(x,y,CLs):
	h, xx, yy = np.histogram2d(x,y,bins=60)
	dx = np.average(np.diff(xx))
	dy = np.average(np.diff(yy))
	h  = gf(h,2).T ; h /= np.matrix.sum(np.matrix(h)*dx*dy)

	h2 = np.sort(h.ravel()) # matrix(n,n) to 1d-array(n*n)
	sums = np.zeros(len(h2))

	for i in range(len(h2)): sums[i] = np.sum(h2[i:])*dx*dy
	func = interp1d(sums,h2,fill_value='extrapolate')

	return xx[1:], yy[1:], h, func(CLs)

def Plot_corner(M, truths, labels, FreeParams, Cov):
	print('\n'+"-----"*10+'\n')
	Np = len(truths)
	print(">> Plotting corner ...")
	plt.rcParams.update({'font.family':'serif'})
	fig = plt.figure(figsize=(2*Np,2*Np))
	for i in range(Np): # columns
		ax = plt.subplot(Np,Np,i*Np+i+1)

		Xi = M[:,i]
		qx = [0.5*(1-0.99),0.5*(1+0.99)]
		li = np.quantile(Xi,qx)

		CL = 0.68
		qs = [0.5*(1-CL),0.5,0.5*(1+CL)]
		q = np.quantile(M[:,i],qs)
		ax.set_title('$%.2f^{+%.2e}_{-%.2e}$' % (q[1],q[1]-q[0], q[2]-q[1]))

		if('d_L' in labels[i]):
			print('-----> DL in params!!!')
			err = 0.5*(q[2]-q[0])
			err *= 100/q[1] 
			ax.plot(np.nan,np.nan,'ko',label="$\Delta d_L$: %.2e" % err + '%')

		deltaY = 0.5*(q[2]-q[0])
		ly1, ly2 = q[1]-3*deltaY, q[1]+3*deltaY

		x, bins = np.histogram(M[:,i],bins=60)
		y = gf(x,2)
		plt.fill_between(bins[1:],np.zeros(len(x)),y,color='k',alpha=0.4)
		plt.plot(bins[1:],y,'k-',lw=2)
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_xlim(ly1,ly2)
		ax.set_ylim(0,1.2*max(y))
		ax.vlines(truths[i],0,1.2*max(y),color='r')
		
		ax.vlines(q[0],0,1.2*max(y),color='k',ls='--')
		ax.vlines(q[1],0,1.2*max(y),color='k',ls='--')
		ax.vlines(q[2],0,1.2*max(y),color='k',ls='--')

		if(i==(Np-1)): ax.set_xlabel(labels[i])
		for j in range(i): # rows
			x = M[:,j]
			y = M[:,i]
			Xj = M[:,j]
			lj = np.quantile(Xj,qs)

			deltaX = 0.5*(lj[2]-lj[0])
			lx1, lx2 = lj[1]-3*deltaX, lj[1]+3*deltaX

			ax = plt.subplot(Np,Np,i*Np+j+1)

			xx, yy, h, levels = Curves_CL(x,y,[0.95,0.68,0.0])
			ax.contour(xx,yy,h,levels=levels,colors='k')
			ax.contourf(xx,yy,h,levels=levels,cmap='Greys')
			
			Cov2 = np.zeros([2,2])
			Cov2[0,0] = Cov[j,j]
			Cov2[1,1] = Cov[i,i]
			Cov2[0,1] , Cov2[1,0] = Cov[j,i], Cov[j,i]

			print("Error(%s): " % labels[j], np.sqrt(Cov2[0][0]) )
			#Add_Ellipse(Cov2,[truths[j],truths[i]],ax=ax, fill=False)

			#-----------------------------------------------------------
			ax.vlines(truths[j],ly1,ly2,color='r')
			ax.hlines(truths[i],lx1,lx2,color='r')
			ax.plot(truths[j],truths[i],'rs',mec='k')

			muY = truths[j]
			sigmaY = np.sqrt(Cov[j,j])
			
			ax.set_xlim(lx1,lx2)
			ax.set_ylim(ly1,ly2)

			ax.set_xticklabels([])
			ax.set_yticklabels([])
			if(i==(Np-1)): 	ax.set_xlabel(labels[j])
			if(j==0): ax.set_ylabel(labels[i])

	ax.plot(np.nan,np.nan,'r--',lw=2,label='From $F^{-1}$')
	print('\n'+"-----"*10+'\n')
	return fig

def gaussian(mu,sigma):
	x1 = mu-3*sigma
	x2 = mu+3*sigma
	x = np.linspace(x1,x2,100)
	y = np.exp(-(x-mu)**2/sigma**2) / np.sqrt(2*np.pi*sigma**2)
	return x, y

def Plot_Fisher(truths, labels, FreeParams, Cov):
	print('\n'+"-----"*10+'\n')
	Np = len(truths)
	print(">> Plotting Fisher ...")
	plt.rcParams.update({'font.family':'serif'})
	fig = plt.figure(figsize=(2*Np,2*Np))
	for i in range(Np): # columns
		ax = plt.subplot(Np,Np,i*Np+i+1)

		mu = truths[i]
		sigma = np.sqrt(Cov[i,i])
		x, y = gaussian(mu,sigma)

		ax.set_title('$%.2f^{+%.2e}_{-%.2e}$' % (mu,sigma,sigma))

		if('DL' in labels[i]):
			try: ax.plot(np.nan,np.nan,'ko',label="$\Delta d_L$ = %.2f" % fisher_err + '%')
			except: pass

		plt.fill_between(x,np.zeros(len(x)),y,color='k',alpha=0.4)
		plt.plot(x,y,'k-',lw=2)
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_xlim(min(x),max(x))
		ax.set_ylim(0,1.2*max(y))
		ax.vlines(truths[i],0,1.2*max(y),color='r')
		
		ax.vlines(mu-sigma,0,1.2*max(y),color='k',ls='--')
		ax.vlines(mu,0,1.2*max(y),color='k',ls='--')
		ax.vlines(mu+sigma,0,1.2*max(y),color='k',ls='--')

		if(i==(Np-1)): ax.set_xlabel(labels[i])

		muX = mu
		sigmaX = sigma
		for j in range(i): # rows
			x, f = gaussian(truths[j],np.sqrt(Cov[j,j]))
			y, f = gaussian(truths[i],np.sqrt(Cov[i,i]))

			ax = plt.subplot(Np,Np,i*Np+j+1)

			Cov2 = np.zeros([2,2])
			Cov2[0,0] = Cov[j,j]
			Cov2[1,1] = Cov[i,i]
			Cov2[0,1] , Cov2[1,0] = Cov[j,i], Cov[j,i]

			Add_Ellipse(Cov2,[truths[j],truths[i]],ax=ax, fill=True)
			
			muY = truths[j]
			sigmaY = np.sqrt(Cov[j,j])

			ax.vlines(truths[j],min(y),max(y),color='r')
			ax.hlines(truths[i],min(x),max(x),color='r')
			ax.plot(truths[j],truths[i],'rs',mec='k')

			ax.set_xlim(muX-3*sigmaX,muX+3*sigmaX)
			ax.set_ylim(muY-3*sigmaY,muY+3*sigmaY)

			ax.set_xticklabels([])
			ax.set_yticklabels([])
			if(i==(Np-1)): 	ax.set_xlabel(labels[j])
			if(j==0): ax.set_ylabel(labels[i])

	print('\n'+"-----"*10+'\n')
	return fig