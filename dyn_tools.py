import numpy as np


def delta_(t,h):
#% Delta Impulse Input
	d=np.copy(np.atleast_1d(t))
	d[d>h] = 0.
	d[d<0.] = 0.
	return d
 
def DxDt_r(x ,t ,p ):
	''''Roessler '''
	if p is False: a,b,c = 0.398,2.,4.0
	else: a,b,c = p[0],p[1],p[2]
	dx = np.zeros(3)
	dx[0] = - x[1] - x[2]
	dx[1] = x[0] + a*x[1]
	dx[2] = b + x[2] *(x[0] - c)
	return dx


def DxDt_mass(x ,t ,p ):
	''''double dumped mass system '''
	w2 = 1.0
	dx = np.zeros(4)
	dx[0] =  x[2] 
	dx[1] = x[3]
	dx[2] = -2*w2*x[0] + w2*x[1]
	dx[3] = w2*x[0] - 2*w2*x[1]
	return dx


def f_ds(x ,t ,p ):
	a,b = -8./7,-5./7
	return b*x + 0.5 * (a - b) * (np.abs(x+1) - np.abs(x-1) )

def DxDt_ds(x ,t ,p ):
	''''Double scroll '''
	if p is False:alpha, beta, gamma = 9.,100./7,0.
	else: alpha, beta, gamma = p[0],p[1],p[2]
	dx = np.zeros(3)
	dx[0] = alpha * ( x[1] - x[0] - f_ds(x[0]) )
	dx[1] = x[0] - x[1] + x[2]
	dx[2] = -beta * x[1] -gamma * x[2]
	return dx

def DxDt_l(x ,t ,p ):
	''''Lorenz '''
	if p is False: sigma, rho, beta = 10., 28., 8./3.
	else: sigma, rho, beta = p[0],p[1],p[2]
	dx = np.zeros(3)
	dx[0] = sigma* (x[1] - x[0])
	dx[1] = x[0]*( rho-x[2])  - x[1]
	dx[2] = x[0] * x[1] - beta* x[2]
	return dx

def DxDt_cord(x ,t ,p ):
	''''Cord attractor'''
	if p is False: a,b,F,G = 0.25,4.0,8.0,1.0	
	else: a,b,F,G = p[0],p[1],p[2], p[3]	
	dx = np.zeros(3)
	dx[0] = -x[1] -x[2] - a*x[0] + a*F
	dx[1] = x[0]*x[1]  - b*x[0]*x[2] - x[1] + G
	dx[2] = b*x[0]*x[1] + x[0]*x[2] - x[2]
	return dx

def DxDt_l84(x ,t ,p ):
	''''Lorenz 84'''
	if p is False: a,b,F,G = 0.25,4.0,8.0,1.0	
	else: a,b,F,G = p[0],p[1],p[2], p[3]	
	dx = np.zeros(3)
	dx[0] = -x[1]*x[1] -x[2]*x[2] - a*x[0] + a*F
	dx[1] = x[0]*x[1]  - b*x[0]*x[2] - x[1] + G
	dx[2] = b*x[0]*x[1] + x[0]*x[2] - x[2]
	return dx


def rk4(f,x,dt,t=False,p=False,u=False):
	''' rk4 integrator
def rk4(f,x,dt,t=False,p=False,u=False)'''

	if t is False: t=0.
	if u is False: ut= lambda t: 0.
	if callable(u) is True: ut = lambda t: u(t)
	else: ut= lambda t: u
	if p is False: fdyn = lambda x,t,p: f(x,t)
	else: fdyn = lambda x,t,p: f(x,t,p)
	k1=fdyn(x,t,p) + ut(t)
	k2=fdyn(x + dt*k1/2,t+dt/2,p) + ut(t+dt/2)
	k3=fdyn(x + dt*k2/2,t+dt/2,p) + ut(t+dt/2)
	k4=fdyn(x + dt*k3,t+dt,p) + ut(t+dt)
	
	return x+dt*(k1+2*k2+2*k3+k4)/6

def integrate(f,x0,dt, t0=0, tfin = 10., p=False,u=False):
	y=[]
	n = int((tfin-t0)/dt)
	time = t0+np.arange(n)*dt
	for it, t in enumerate(time):
		if it==0:y.append(x0)
		else:y.append(rk4(f,y[-1],dt,t,p,u))
	return y,time.tolist()
