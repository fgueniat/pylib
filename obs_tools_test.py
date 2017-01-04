
import numpy as np
import scipy
from la_tools import sparse
import dyn_tools as dyn
# if available import pylab (from matlibplot)
try:
	import matplotlib.pylab as plt
except ImportError:
	pass

def emgr(f, g, s, t, w, pr = 0., nf = 0., ut = 1., us = 0., xs = 0., um = 1., xm = 1., pm =np.array([]), DOT = False, ODE = False):

	# Local Variables: a0, tu, xm, WM, up, d, xx, WT, o, xs, DOT, pr, pp, TX, TU, nf, K, ODE, pm, A, C, D, G, F, i1, i0, DWX, M, L, N, Q, P, n0, W, V, ux, c, tx, g, f, uu, ut, k, m, us, mh, n, um, p, s, sh, u, t, w, h, y, x, z
	# Function calls: scales, emgr, avg, find, size, ainv, floor, sum, sqrt, any, zeros, pscales, pi, max, cos, ceil, nargin, ones, isempty, numel, isa, lower, min, sparse, error, Inf, strcmp
	#%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
	#%%% version: 5.0 ( 20150-20 )
	#%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
	#%%% license: BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
	#%
	#%% SYNTAX:
	#%   W = emgr(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[pm]);
	#%
	#%% SUMMARY:
	#%   empirical gramian matrix computation for model reduction,
	#%   decentralized control, sensitivity analysis, parameter identification,
	#%   uncertainty quantification and combined state and parameter reduction
	#%   of large-scale input-output systems. Enables data-driven analysis of
	#%   input-output coherence and gramian-based nonlinear model order reduction.
	#%   Compatible with OCTAVE and MATLAB.
	#%
	#%% ARGUMENTS:
	#%   (handle) f : system function handle; signature: xdot = f(x,u,p,t)
	#%   (handle) g : output function handle; signature:	y = g(x,u,p,t)
	#%   (vector) s : system dimensions [inputs,states,outputs]
	#%   (vector) t : time discretization [step,stop]
	#%   (string) w : character encoding gramian type:
	#%	  * 'c' empirical controllability gramian (WC)
	#%	  * 'o' empirical observability gramian (WO)
	#%	  * 'x' empirical cross gramian (WX aka WCO and XCG)
	#%	  * 'y' empirical linear cross gramian (WY)
	#%	  * 's' empirical sensitivity gramian (WS)
	#%	  * 'i' empirical identifiability gramian (WI)
	#%	  * 'j' empirical joint gramian (WJ)
	#%   (matrix) pr = 0 : parameters, each column is one set
	#%   (vector) nf = 0 : option flags, ten components:
	#%	  * centering: no(0), init(1), steady(2), mean(3), rms(4), midrange(5)
	#%	  * input scales: single(0), linear(1), geometric(2), log(3), sparse(4)
	#%	  * state scales: single(0), linear(1), geometric(2), log(3), sparse(4)
	#%	  * input rotations: unit(0), single(1)
	#%	  * state rotations: unit(0), single(1)
	#%	  * scale type: no(0), Jacobi preconditioner(1), steady-state(2)
	#%	  * cross gramian type (only: WX,WY,WJ): regular(0), non-symmetric(1) (WZ)
	#%	  * extra input (only: WO,WX,WS,WI,WJ): no(0), param(1), state(2), both(3)
	#%	  * parameter centering (only: WS,WI,WJ): no(0), linear(1), logarithmic(2)
	#%	  * Schur-complement (only: WI,WJ): detailed(0), approximate(1)
	#%   (handle) ut = 1 : input function handle; default: delta impulse
	#%   (vector) us = 0 : steady-state input
	#%   (vector) xs = 0 : steady-state and initial state x0
	#%   (matrix) um = 1 : input scales
	#%   (matrix) xm = 1 : initial-state scales
	#%   (matrix) pm = [] : parameter scales (reserved)
	#%
	#%% RETURNS:
	#%   (matrix) W : Gramian Matrix (only: WC,WO,WX,WY)
	#%	 (cell) W : [State-, Parameter-] Gramian (only: WS,WI,WJ)
	#%
	#%% CITATION:
	#%   C. Himpe (2016). emgr - Empirical Gramian Framework (Version 5.0)
	#%   [Software]. Available from http://gramian.de . doi: 10.5281/zenodo.162135 .
	#%
	#%% SEE ALSO:
	#%   gram
	#%
	#%% KEYWORDS:
	#%   model reduction, empirical gramians, cross gramian matrix, MOR
	#%
	#% Further information: <http://gramian.de>
	#%$
	#% Inner Product Handle
	#% Integrator Handle
	#% Distributed Cross Gramian (Column Width, Index)
	if DOT is False:DOT = np.dot
#
	if ODE is False:ODE = ssp2
#
	#% Version Info
	if f == 'version':
		W = 5.0
		return W
#
########################################################################################################################
########################################################################################################################
#					#%% GENERAL SETUP
########################################################################################################################
########################################################################################################################
#
	#% System Dimensions
	M = s[0] #% Number of inputs
	N = s[1] #% Number of states
	Q = s[2] #% Number of outputs
	if len(s) == 4.: #% Number of augmented parameter-states / inputs
		A = s[3]
	else:
		A = 0.
#
	pr = np.atleast_2d(pr)
	P = pr.shape[0] #% Dimension of parameter (each column is a set of parameter)
	K = pr.shape[1] #% Number of parameters
	h = t[0] #% Width of time step
	L = int(t[1]/h) + 1. #% Number of time steps plus initial value

	w = w.lower() #% Ensure lower case gramian type

	#% Lazy Arguments
	if g == 1. or g is False: # observable function is identity function
		g = lambda x, u, p, t: x
		Q = N
#
	print(nf)
	nf = np.atleast_1d(nf)
	if len(nf)<10.:
		nf = np.r_[nf,np.zeros(10-len(nf))]
#		

#############################################
#############################################
# to be fixed:
# to be fixed:
#############################################
#############################################
	ut = np.atleast_1d(ut)
	if len(ut) == 1.:
		#% and shoudl be an array ?
		if ut == np.inf:
		#% Linear Chirp Input
			mh = np.ones((M, 1.))
			sh = (1.0/L - 0.1 / h) / L
			ut = lambda t: 0.5 + np.dot( 0.5 * mh * np.cos(np.dot(np.dot(2.0, np.pi), (matdiv(0.1, h)+np.dot(np.dot(0.5, sh), t))*t)))
		else:
		#% Delta Impulse Input
			mh = np.ones((M, 1.))/h
			ut = lambda t: mh * dyn.delta_(t,h)
	elif callable(ut) is False:
	#% Delta Impulse Input
		mh = np.ones((M, 1.))/h
		ut = lambda t: mh * dyn.delta_(t,h)

#############################################
#############################################
#############################################
#############################################
#	
	us = np.atleast_1d(us)
	if len(us) == 1.: us = us * np.ones((M, 1.))
#
	xs = np.atleast_1d(xs)
	if len(xs) == 1.: xs = xs * np.ones((N, 1.))
#
	um = np.atleast_2d(um)
	if len(um) == 1.: um = np.atleast_2d(um * np.ones((M, 1.)))
#
	xm = np.atleast_1d(xm)
	if len(xm) == 1.: 
		if w == 'y':
			xm =   xm * np.ones((M, 1.))
		else:
			xm = xm * np.ones((N, 1.))
#
	print(um)
	if um.shape[1] == 1.: um = scales(um, nf[1], nf[3])
#
	if xm.shape[1] == 1.: xm = scales(xm, nf[2], nf[4])
#
	#%% GRAMIAN SETUP
	if w == 'c' or w == 'o' or w == 'x' or w == 'y':
		C = um.shape[1]
		#% Number of input scales
		D = xm.shape[1]
		#% Number of state scales
		if isempty(pm): pm = np.zeros((1., np.max([C, D])))
#
		switch_val=nf[5]
		if False: # switch 
			pass
		elif switch_val == 1.:
			#% Preconditioned run
			nf[5] = 0.
			WT = emgr(f, g, s, t, w, pr, nf, ut, us, xs, um, xm, pm)
			TX = np.sqrt(WT[0::N+1.]).conj().T
			tx = 1.0/TX
			F = f
			f = lambda x, u, p, t: TX*F(np.dot(tx,x) ,u ,p ,t )
			G = g
			g = lambda x, u, p, t: G(np.dot(tx,x) ,u ,p ,t )
		elif switch_val == 2.:
			#% Steady-state scaled run
			TU = us[:,0]
			TU[TU == 0.] = 1.0
			tu = 1.0/TU
			TX = xs
			TX[TX == 0.] = 1.0
			tx = 1.0/TX
			F = f
			f = lambda x, u, p, t: TX*F(np.dot(tx,x),np.dot(tu,u) ,p ,t)
			G = g
			g = lambda x, u, p, t: G(np.dot(tx,x) ,np.dot(tu,u) ,p ,t )
		
		if nf[7] == 1. or nf[7] == 3.:
			#% Extra input for parameter perturbations
			up = lambda t: us+ut(t)
		else:
			up = lambda t: us
			
		
		if nf[7] == 2. or nf[7] == 3.:
			#% Extra input for state perturbations
			ux = lambda t: us+ut(t)
		else:
			ux = lambda t: us
			
		
	
	
	#%% GRAMIAN COMPUTATION
	switch_val=w
	if False: # switch 
		pass
	#% Empirical gramian types
	elif switch_val == 'c':
		#% Controllability gramian
		W = np.zeros((N, N))
		o = np.zeros((A, 1.))
		for k in np.arange(1., K+1):
			for c in np.arange(1., C+1):
				for m in nonzero(um[:,c]).conj().T:
					uu = lambda t: us + np.dot(ut(t),sparse(m,1,um(m,c),M,1))
					####################################
					# to fix
					print('to fix')
					####################################
					x = ODE(f,False,t,xs,uu,pr[:,k])
					####################################
					####################################
					x = x - avg(x,nf[0])
					x = np.dot( x, (1.0/um(m,c)) )
					W = W + np.dot(x,x.conj().T)
				for m in nonzero(um[:,c]).conj().T:
					pp = pr[:,k] + sparse(m,1,pm(m,c),P,1);
					####################################
					# to fix
					print('to fix')
					####################################
					x = ODE(f,False,t,xs,up,pp);
					x = x - avg(x,nf[1])
					x = np.dot( sx, (1.0/pm(m,c)) )
					o[m] = o[m] + np.sum(np.dot(x,x))
#
		W = W* h/( np.dot(C, K) ) 
		W = 0.5* (W+W.conj().T)
		print('to fix 221 ')
		if A > 0.: W = [W,o/np.sum(o)]


# OSEF les autres gram
#	elif switch_val == 'o':
#		#% Observability gramian
#		W = np.zeros((N+A,N+A))
#		o = np.zeros((Q*L, N+A))
#		for k in np.arange(1., K+1):
#			for d in np.arange(1., D+1):
#				for n in nonzero(xm[:,d]).conj().T:
#					#% parfor
#					xx = xs + sparse(n,1,xm(n,d),N,1)
#					y = dyn.integrate(???)
#		#                       y = ODE(f,g,t,xx,ux,pr(:,k));
#					y = y - avg(y,nf[1])
#					y = y  * (1.0/xm(n,d))
#					o[:,n] = y[:]
############## check point
#				for n in nonzero(pm[:,d]).conj().T:
#					#% parfor
#					pp = pr[:,k] + sparse(n,1,pm[n,d],P,1);
#					y = np.integrate(f??);
##					y = ODE(f,g,t,xs,up,pp);
#					y = y - avg(y,nf[1]);
#					y = y * (1.0/pm[n,d]);
#					o[:,N+n] = y[:];
##
#				W = W + np.dot(o.conj().T,o)
##
#		W = W * (h/(D*K))
#		W = np.dot(0.5, W+W.conj().T)

#	elif switch_val == 'x':
#		#% Cross gramian
#		if M != Q:
#			if nf[6] == 0.:
#				matcompat.error('emgr: non-square system!')
#			
#			
#		
#		
#		if isempty(DWX):
#			#% Full cross gramian
#		n0 = 0.
#		a0 = -N
#		W = np.zeros(N, (N+A))
#		o = np.zeros(L, (N+A), Q)
#		else:
#			#% Distributed cross gramian
#			if np.any((np.ceil(DWX) != np.floor(DWX))):
#				matcompat.error('emgr: non-positive or non-integer values in DWX!')
#			
#			
#			if np.any((DWX<=0.)):
#				matcompat.error('emgr: non-positive or non-integer values in DWX!')
#			
#			
#			i0 = np.dot(DWX[1]-1., DWX[0])+1.
#			if i0<=N:
#				#% State-space columns setup
#			i1 = matcompat.max((i0+DWX[0]-1.), N)
#			xm[np.array(np.hstack((0:i-1., int(i1+1.)-1:))),:] = 0.
#			pm = np.zeros(1., D)
#			n0 = i-1.
#			else:
#				#% Parameter-space columns setup
#				i0 = i0-np.dot(np.ceil(matdiv(N, DWX[0])), DWX[0])-N
#				i1 = matcompat.max((i0+DWX[0]-1.), (N+A))
#				xm = np.zeros(1., D)
#				pm[np.array(np.hstack((0:i0-N-1., int(i1-N+1.)-1:))),:] = 0.
#				a0 = i0-N-1.
#				
#			
#			W = np.zeros(N, (i1-i0+1.))
#			o = np.zeros(L, (i1-i0+1.), Q)
#			if i0 > i1:
#				W = 0.
#				return []
#			
#			
#			
#		
#		for k in np.arange(1., (K)+1):
#			for d in np.arange(1., (D)+1):
#				for n in nonzero(xm[:,int(d)-1]).conj().T:
#					#% parfor
#					
#				for n in nonzero(pm[:,int(d)-1]).conj().T:
#					#% parfor
#					
#				if nf[6]:
#					#% Non-symmetric cross gramian: cache average
#				o[:,:,0] = np.sum(o, 3.)
#				
#				for c in np.arange(1., (C)+1):
#					for m in nonzero(um[:,int(c)-1]).conj().T:
#						#% parfor
#						
#					
#				
#			
#		W = np.dot(W, matdiv(h, np.dot(np.dot(C, D), K)))
#	elif switch_val == 'y':
#		#% Linear cross gramian
#		if M != Q:
#			if nf[6] == 0.:
#				matcompat.error('emgr: non-square system!')
#			
#			
#		
#		
#		W = np.zeros(N, N)
#		o = np.zeros(N, L, M)
#		for k in np.arange(1., (K)+1):
#			for c in np.arange(1., (C)+1):
#				for m in nonzero(um[:,int(c)-1]).conj().T:
#					#% parfor
#					
#				if nf[6]:
#					#% Non-symmetric cross gramian: cache average
#				o[:,:,0] = np.sum(o, 3.)
#				
#				for m in nonzero(xm[:,int(c)-1]).conj().T:
#					#% parfor
#					
#				
#			
#		W = np.dot(W, matdiv(h, np.dot(C, K)))
#	elif switch_val == 's':
#		#% Sensitivity gramian
#		[pr, pm] = pscales(pr, matcompat.size(um, 2.), nf[8])
#		W = emgr(f, g, np.array(np.hstack((M, N, Q, P))), t, 'c', pr, nf, ut, us, xs, um, xm, pm)
#		#% W{1} % Controllability gramian
#		#% W{2} % Sensitivity gramian diagonal
#	elif switch_val == 'i':
#		#% Identifiability gramian
#		[pr, pm] = pscales(pr, matcompat.size(xm, 2.), nf[8])
#		V = emgr(f, g, np.array(np.hstack((M, N, Q, P))), t, 'o', pr, nf, ut, us, xs, um, xm, pm)
#		W.cell[0] = V[0:N,0:N]
#		#% Observability gramian
#		WM = V[0:N,int(N+1.)-1:N+P]
#		if nf[9]:
#			#% Identifiability gramian
#		W.cell[1] = V[int(N+1.)-1:N+P,int(N+1.)-1:N+P]
#		else:
#			W.cell[1] = V[int(N+1.)-1:N+P,int(N+1.)-1:N+P]-np.dot(np.dot(WM.conj().T, ainv(W.cell[0])), WM)
#			
#		
#	elif switch_val == 'j':
#		#% Joint gramian
#		[pr, pm] = pscales(pr, matcompat.size(xm, 2.), nf[8])
#		V = emgr(f, g, np.array(np.hstack((M, N, Q, P))), t, 'x', pr, nf, ut, us, xs, um, xm, pm)
#		if isempty(DWX) == 0.:
#			#% Distributed cross gramian
#		W = V
#		return []
#		
#		W.cell[0] = V[0:N,0:N]
#		#% Cross gramian
#		WM = V[0:N,int(N+1.)-1:N+P]
#		if nf[9]:
#			#% Cross-identifiability gramian
#		W.cell[1] = np.dot(-0.5, np.dot(WM.conj().T, WM))
#		else:
#			W.cell[1] = np.dot(-0.5, np.dot(np.dot(WM.conj().T, ainv((W.cell[0]+W.cell[0].conj().T))), WM))
#			
#		
	
	else:
		print('emgr: unknown gramian type!')
		W = -1
	#%% ======== INPUT AND STATE SCALES ========
	return [W]





def scales(s, d, e):

	# Local Variables: s, e, d
	# Function calls: scales
	#%%% summary: scales (Input and initial state perturbation scales)
	#%$
	switch_val=d
	if False: # switch 
		pass
	elif switch_val == 1.:
		#% Linear
		s = np.dot(s, np.array(np.hstack((0.25, 0.50, 0.75, 1.0))))
	elif switch_val == 2.:
		#% Geometric
		s = np.dot(s, np.array(np.hstack((0.125, 0.25, 0.5, 1.0))))
	elif switch_val == 3.:
		#% Logarithmic
		s = np.dot(s, np.array(np.hstack((0.001, 0.01, 0.1, 1.0))))
	elif switch_val == 4.:
		#% Sparse
		s = np.dot(s, np.array(np.hstack((0.38, 0.71, 0.92, 1.0))))
	
	if e == 0.:
		s = np.array(np.hstack((-s, s)))
	
	
	#%% ======== PARAMETER SCALES ========
	return [s]
def pscales(p, n, e):

	# Local Variables: pr, e, n, pmin, p, lmax, lmin, pmax, pm
	# Function calls: real, log, min, max, linspace, pscales, error, exp, size
	#%%% summary: pscales (Parameter perturbation scales)
	#%$
	if matcompat.size(p, 2.) == 1.:
		matcompat.error('emgr: min + max parameter required!')
	
	
	pmin = matcompat.max(p, np.array([]), 2.)
	pmax = matcompat.max(p, np.array([]), 2.)
	switch_val=e
	if False: # switch 
		pass
	#% Parameter centering
	elif switch_val == 1.:
		#% Linear
		pr = np.dot(0.5, pmax+pmin)
		pm = np.dot(pmax-pmin, np.linspace(0., 1.0, n))
		pm = pm+pmin-pr
	elif switch_val == 2.:
		#% Logarithmic
		lmin = np.log(pmin)
		lmax = np.log(pmax)
		pr = np.real(np.exp(np.dot(0.5, lmax+lmin)))
		pm = np.dot(lmax-lmin, np.linspace(0., 1.0, n))
		pm = pm+lmin
		pm = np.real(np.exp(pm))+pmin-pr
	else:
		#% None
		pr = pmin
		pm = np.dot(pmax-pmin, np.linspace(0., 1.0, n))
	
	#%% ======== TRAJECTORY AVERAGE ========
	return [pr, pm]
def avg(x, d):

	# Local Variables: x, m, d
	# Function calls: size, min, max, sum, sqrt, zeros, avg, mean
	#%%% summary: avg (State and output trajectory centering)
	#%$
	switch_val=d
	if False: # switch 
		pass
	elif switch_val == 1.:
		#% Initial state / output
		m = x[:,0]
	elif switch_val == 2.:
		#% Steady state / output
		m = x[:,int(0)-1]
	elif switch_val == 3.:
		#% Mean state / output
		m = np.mean(x, 2.)
	elif switch_val == 4.:
		#% Root-mean-square state / output
		m = np.sqrt(np.sum((x*x), 2.))
	elif switch_val == 5.:
		#% Midrange state / output
		m = np.dot(0.5, matcompat.max(x, np.array([]), 2.)-matcompat.max(x, np.array([]), 2.))
	else:
		#% None
		m = np.zeros(matcompat.size(x, 1.), 1.)
	
	#%% ======== FAST APPROXIMATE INVERSION ========
	return m
def ainv(m):

	# Local Variables: x, m, d, n
	# Function calls: diag, ainv, numel
	#%%% summary: ainv (Quadratic complexity approximate inverse matrix)
	#%$
	d = np.diag(m)
	d[int((d != 0.))-1] = 1.0/d[int((d != 0.))-1]
	n = len(d)
	x = m*-d
	x = x*d.conj().T
	x[0::n+1.] = d
	#%% ======== DEFAULT ODE INTEGRATOR ========
	return x


def ssp2(f,g,t,x0,u=False,p=False):
	''' integrator:
		f dynamical function: dxdt=f(x,t) + u(t)
		u control function
		g observable function
		p parameters
	'''
		
	if callable(g) is not True:
		g = lambda x,ut,p,t: x
	if isinstance(u,(float,int)) is True:
		ut = lambda t: u
	if callable(u) is False or u is False:
		ut = lambda t: 0
	else: ut = lambda t: u(t)
	dt = t[0]
	N = int(t[1]/dt) + 1

	x = dyn.integrate(f,x0,dt,p=p,tfin = N*dt, u=ut)
	time,x = x[1],x[0]
	y = [g(x[i],ut(time[i]),p,time[i]) for i in xrange(len(x))]
	return y

def ssp2_bis(f, g, t, x0, u, p):

	# Local Variables: xk, g, f, h, K, p, s, u, t, uk, hs, y, x, x0, STAGES, k, tk
	# Function calls: isscalar, ssp2, floor
	#%%% summary: ssp2 (Low-Storage Stability Preserving Runge-Kutta SSP32)
	#%$
	if np.isscalar(STAGES) == 0.:
		STAGES = 3.
	
	
	h = t[0]
	K = np.floor(matdiv(t[1], h))+1.
	x = x0
	y = g[int(x)-1,int(u[-1])-1,int(p)-1,-1]
	y[int(0)-1,int(K)-1] = 0.
	#% Preallocate trajectory
	hs = matdiv(h, STAGES-1.)
	xk = x
	for k in np.arange(1., (K-1.)+1):
		tk = np.dot(k, h)
		uk = u[int(tk)-1]
		for s in np.arange(1., (STAGES-1.)+1):
			xk = xk+np.dot(hs, f[int(xk)-1,int(uk)-1,int(p)-1,int(tk)-1])
			tk = tk+hs
			
		xk = (np.dot(STAGES-1., xk)+x+np.dot(h, f[int(xk)-1,int(uk)-1,int(p)-1,int(tk)-1]))/STAGES
		x = xk
		y[:,int((k+1.))-1] = g[int(x)-1,int(uk)-1,int(p)-1,int(tk)-1]
		
	return y
