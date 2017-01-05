import numpy as np
from scipy.optimize import minimize as solve
from scipy.linalg import sqrtm
import dyn_tools as dyn
import plot_tools as pt
from cantera_tools import chemreac

measures_possible = ['cond', 'det','trace']

def slow_manifold(alpha_0, dx=False, f_def_obs=False,eps_ = False,time = False):
	if dx is False:
		print 'run f_explore first'
		return -1
	if time is False:
		print 'time needed'
		return -1
	if eps_ is False:
		print 'need precision'
		return -1
	if f_def_obs is False:
		print ' obserabl def as y=sum alpha_i x_i'
		f_def_obs = lambda x,alpha: np.dot(alpha,x)
	alpha = solve(cost_gramian,alpha_0,args = (f_def_obs,0.,dx,eps_,time))
	return alpha

def cost_gramian(alpha, f_def_obs=False, rho=0., dx=False, eps_=False, time=False,nt = False, ny=False, nx = False):
	#def obs
	if f_def_obs is False:fobs = lambda xx: xx[1] - np.sum([xx[0]**a for a in alpha])
	else:fobs = lambda x: f_def_obs(x,alpha)
	if nt is False:nt = np.size(time)
	if nx is False:nx = np.size(dx[0,:,0,0])
	if ny is False:ny = np.size(fobs(dx[0,:,0,0]))
	w = EG(fobs = fobs, dx=dx, eps_ = eps_, time = time, nt = nt, nx = nx, ny = ny)
	J =1./ M(w) + rho * np.linalg.norm(alpha)
	return J

def EG(fobs,dx,eps_,time,nt = False,nx = False,ny = False,offset = 1./2):
	''' Empirical gramian'''
	if nt is False:nt = np.size(time)
	if nx is False:nx = dx[0,:,0,0].size
	if ny is False:ny = fobs(dx[0,:,0,0]).size
	offset = int(offset*nt)+1
	W = np.zeros((nx,nx))
	for it in xrange(offset,nt):
		dt = time[it]-time[it-1]
		y=np.zeros((nx,2*ny))
		for ix in xrange(nx):
			y[ix,0:ny] = fobs(dx[it,:,ix,0]);
			y[ix,ny:] = fobs(dx[it,:,ix,1]);
		W = W+dt*np.dot(y,y.T) # equiv to W = W + y*y.T

	W = W/(4.*eps_*eps_);
	W = W+W.T
	return W

def ode(fdyn,time,x0):
	dt = time[1]-time[0]
	nt = time.size
	nx = x0.size
	X = np.zeros((nt,nx))
	for it,t in enumerate(time):
		if it==0: X[it,:] = x0
		else: X[it,:] = dyn.rk4(fdyn,X[it-1,:],dt,t)
	return X

def f_explore(fdyn,x0,time,eps_ = False):
	'''compute the perturbation around a given trajectory  defined by the
initial conditions x0. 
output has dimensions (nt x nx) x (nx x 2).
(nt x nx) trajectories
(nx x 2) positive and negative perturbation around each component.'''
	
	nx = x0.size
	nt = time.size
	
	dx = np.zeros((nt,nx,nx,2))
	if isinstance(fdyn,chemreac):
		for i in xrange(nx):
			base = np.zeros(nx)
			base[i] = eps_
			reaction_p = chemreac(cti=fdyn.cti,V=fdyn.V,problem = fdyn.problem,concentrations=x0+base)
			reaction_p.integrate(time=time)
			xp = np.array(reaction_p.h_z)[1:,:]
			reaction_m = chemreac(cti=fdyn.cti,V=fdyn.V,problem = fdyn.problem,concentrations=x0-base)
			reaction_m.integrate(time=time)
			xm = np.array(reaction_m.h_z)[1:,:]
			dx[:,:,i,0] = xp
			dx[:,:,i,1] = xm


	else:
		#x = ode(fdyn,time,x0)

		for i in xrange(nx):
			base = np.zeros(nx)
			base[i] = eps_
			xp = ode(fdyn,time,x0+base)
			xm = ode(fdyn,time,x0-base)
			dx[:,:,i,0] = xp
			dx[:,:,i,1] = xm

	return dx

def p_explore(dd):
	''' plot the results of f_explore'''
	nt = dd.shape[0]
	nx = dd.shape[1]
	xx = ()
	yy = ()
	for i in xrange(nx):
		xx=xx+(dd[:,0,i,0],)
		yy=yy+(dd[:,1,i,0],)
		xx=xx+(dd[:,0,i,1],)
		yy=yy+(dd[:,1,i,1],)
	pd=pt.Paradraw()
	pd.x_label = '$x_1$'
	pd.y_label = '$x_2$'
	pd.marks = ['k']	
	pt.multiplot2(xx,yy,pd)
	return xx,yy


def M(W,measure='cond'):
	''' Measure of the quality of the Gramian '''

	if measure is 'cond':
		J = np.linalg.cond(W)
	if measure is 'det':
		try:
			J=np.linalg.det( np.linalg.inv(W) )
		except:
			J = 1.e18
	if measure is 'trace':
		try:
			J=np.trace( np.linalg.inv(W) )
		except:
			J = 1.e18
	if J == np.inf:
		J = 1.e18
	return np.log(J)


def cost_neighboor(alpha,Wi, alpha_0):	

	W = Wi[0]
	for i in xrange(len(alpha)): W += alpha[i]*Wi[i+1]
	J = - M(W,'cond') + 0.*np.linalg.norm(alpha-alpha_0)
	# check tangent ?
	return J

	dx0 = np.zeros((2,n,n))
