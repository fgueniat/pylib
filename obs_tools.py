import numpy as np
from scipy.optimize import minimize as solve
from scipy.linalg import sqrtm
import dyn_tools as dyn
import plot_tools as pt
from cantera_tools import chemreac

measures_possible = ['cond', 'det','trace']

def point2manifold(alpha,x):
	pass

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
	dt = time[1]-time[0]
	offset = int(offset*nt)
	W = np.zeros((nx,nx))
	for it in xrange(offset,nt):
		y=np.zeros((nx,2*ny))
		for ix in xrange(nx):
			y[ix,0:ny] = fobs(dx[it,:,ix,0]);
			y[ix,ny:] = fobs(dx[it,:,ix,1]);
		W = W+np.dot(y,y.T) # equiv to W = W + y*y.T

	W = dt*W/(4.*eps_*eps_);
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


def emp_gram(DY,eps):
	if len(DY[0].shape)==1:DY = [dy[np.newaxis] for dy in DY]
	W = np.array(1./(4*eps**2) * sum([np.dot(dy.T,dy) for dy in DY]))
#	return W
	return sqrtm(W,disp=False)[0]

def obs(C,x):
	return np.dot(C,x)

def obs_matrix(i,j,n,m):
	''' Build the observability matrix '''
	C = np.zeros((n,m))
	for ind_,i_i in enumerate(i):
		C[i[ind_],j[ind_]] =1
	return C.T

def M(W,measure='trace'):
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
		xx=xx+(dd[:,0,i,1],)
		yy=yy+(dd[:,1,i,1],)
	pd=pt.Paradraw()
	pd.x_label = '$x_1$'
	pd.y_label = '$x_2$'
	pd.marks = ['k']	
	pt.multiplot2(xx,yy,pd)
	return xx,yy


def emp_gram(DY,eps):
	if len(DY[0].shape)==1:DY = [dy[np.newaxis] for dy in DY]
	W = np.array(1./(4*eps**2) * sum([np.dot(dy.T,dy) for dy in DY]))
#	return W
	return sqrtm(W,disp=False)[0]

def obs(C,x):
	return np.dot(C,x)

def obs_matrix(i,j,n,m):
	''' Build the observability matrix '''
	C = np.zeros((n,m))
	for ind_,i_i in enumerate(i):
		C[i[ind_],j[ind_]] =1
	return C.T

def M(W,measure='trace'):
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
	return J

def cost(alpha,Wi, measures = ['trace']):
	''' Cost function '''
	rho = 0.
	W = Wi[0]
	for i in xrange(len(alpha)): W += alpha[i]*Wi[i+1]
	J = rho*np.linalg.norm(alpha,1)
	for measure in measures_possible:
		if measure in measures:
			J += M(W,measure=measure)
#	if J<0:J=1.e6
#	print(np.abs(J))
	return J

def cost2(alpha,Wi, measures = ['cond']):
	''' Cost function '''
	rho = 0.

	for i in xrange(len(alpha)):
		if i==0: W = alpha[i]*Wi[i]
		else: W += alpha[i]*Wi[i]
	J = rho*np.linalg.norm(alpha,1)
	for measure in measures_possible:
		if measure in measures:
			J += M(W,measure=measure)
#	if J<0:J=1.e6
#	print(np.abs(J))
	return J

def cost_neighboor(alpha,Wi, alpha_0):	

	W = Wi[0]
	for i in xrange(len(alpha)): W += alpha[i]*Wi[i+1]
	J = - M(W,'cond') + 0.*np.linalg.norm(alpha-alpha_0)
	# check tangent ?
	return J

def cost_neighboor2(alpha,Wi, alpha_0):	
	for i in xrange(len(alpha)):
		if i==0: W = alpha[i]*Wi[i]
		else: W += alpha[i]*Wi[i]
	J = - M(W,'cond') + 0.*np.linalg.norm(alpha-alpha_0)
	# check tangent ?
	return J


def delta_obs(trajs,C):
	n = trajs.shape[1]
	N = trajs.shape[3]
	DY = []
	for it in xrange(N):
		dy = np.array( [  obs(C,trajs[0,ix,:,it]) for ix in xrange(n)] ).flatten() - np.array( [obs(C,trajs[1,ix,:,it]) for ix in xrange(n)] ).flatten()
		DY.append(dy)
	return DY

def opt_gram(func=False,x_0 = False, eps = False,t_0=0,dt=1e-6,tmax = False, bds = False):
# construction of the data:

	if func is False or x_0 is False:
		print "need parameters"
		return -1
	if eps is False: eps = 1.e-3*np.linalg.norm(x_0)
	if tmax is False: tmax = 1000*dt
	if t_0>tmax:
		print "t0/tmax incorrect"
		return -1

	n = len(x_0)
	if len(bds)!=n and len(bds)!=1 and bds is not False:
		print len(bds)
		print "bds not well defined"
		return -1
	N = int((tmax-t_0)/dt)

	time = np.arange(N)*dt
	trajs = np.zeros((2,n,n,N))
	dx0 = np.zeros((2,n,n))
