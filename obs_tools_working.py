import numpy as np
from scipy.optimize import minimize as solve
from scipy.optimize import basinhopping as recuit
from scipy.linalg import sqrtm
import dyn_tools as dyn
import plot_tools as pt
from cantera_tools import chemreac

measures_possible = ['cond', 'det','trace']

def slow_manifold(alpha_0, dx=None, f_def_obs=None,eps_ = None,time = None,rho = 1.e-3, measure = None, verbose = False,offset = 0.2,method = 'BFGS',isrecuit=True, niter=50,maxiter=[2000,250],n_recuit=50):
        '''Identify geometric constra ins based on observability analysis'''
        if dx is None:
                print('run f_explore first')
                return -1
        if time is None:
                print('time needed')
                return -1
        if eps_ is None:
                print('need precision')
                return -1
        if f_def_obs is None:
                print(' obserabl def as y=sum alpha_i x_i')
                f_def_obs = lambda x,alpha: np.dot(alpha,x)
        argums = (f_def_obs,rho,dx,eps_,time,False,False,False,measure,offset,verbose)
        if type(maxiter) is not list: 
                if type(maxiter) is int or type(maxiter) is float:
                        maxiter=[maxiter]
                else:
                        print('maxiter ill defined')
                        return -1
        if isrecuit is True:#Annealing
                opt = {'maxiter':maxiter[1]}
                min_kwargs = {'args':argums,'method':method,'options':opt}
                np.random.seed(111111)
                alpha_0 = recuit(cost_gramian,alpha_0,niter=n_recuit,minimizer_kwargs = min_kwargs,T=10.)
                np.random.seed()
                alpha_0 = alpha_0.x
        opt = {'maxiter':maxiter[0]}
        #if no annealing, we finish with a regular minimization
        a = solve(cost_gramian,alpha_0,method = method, tol = 1.e-6,args = argums,options=opt)
        return a

def cost_gramian(alpha, f_def_obs=False, rho=0., dx=False, eps_=False, time=False, nt = False, ny=False, nx = False, measure=None, offset = 0.2, verbose = False):
        '''Compute the gramian and its  measure'''
        # parameters
        if f_def_obs is False:#def obs
                fobs = lambda xx: xx[1] - np.sum([xx[0]**a for a in alpha])
        else:fobs = lambda x: f_def_obs(x,alpha)
        if nt is False:nt = np.size(time)
        if nx is False:nx = np.size(dx[0,:,0,0])
        if ny is False:ny = np.size(fobs(dx[0,:,0,0]))
        #
        #computation of the gramian
        w = EG(fobs = fobs, dx=dx, eps_ = eps_, time = time, nt = nt, nx = nx, ny = ny, offset = offset)
        if measure is None or measure not in measures_possible:
                measure = 'trace'
        norm_type = 2   
        if verbose is True:
                print('gramian')
                print(w)
                print('Measure of the gramian:' + str(1./(M(w,measure=measure))))
                print('norm: ' + str(np.linalg.norm(alpha,ord=norm_type)) + ' and rho: ' + str(rho))
                print('regularization: ' + str(np.linalg.norm(alpha,ord=norm_type) * rho))

        #
        # cost
        J =1./ M(w,measure = measure) + rho * np.linalg.norm(alpha,ord=norm_type)
        return J


def EG(fobs,dx,eps_,time,nt = False,nx = False,ny = False,offset = .2):
        ''' Empirical gramian'''
        # parameters
        if nt is False:nt = np.size(time)
        if nx is False:nx = dx[0,:,0,0].size
        if ny is False:ny = fobs(dx[0,:,0,0]).size
        offset = int(offset*nt)+1
        offtype = 'slow'
        if offtype == 'fast':
                nstart = 1
                nend = offset
        elif offtype == 'slow':
                nstart = offset
                nend = nt
        else:
                nstart = 1
                nend = nt

        #initialization
        W = np.zeros((nx,nx))
        #main loop for the empirical gramian
        for it in range(nstart,nend):
                dt = time[it]-time[it-1]
                y=np.zeros((nx,2*ny))
                for ix in range(nx):# observation
                        if type(eps_) is float:epsi=2*eps_
                        else:epsi=2*eps_[ix]
                        y[ix,0:ny] = fobs(dx[it,:,ix,0]).flatten()/epsi
                        y[ix,ny:] = fobs(dx[it,:,ix,1]).flatten()/epsi
                # contribution of each observation to the gramian       
                W = W+dt*np.dot(y,y.T) # equiv to W = W + y*y.T
        #re normalization
        #W = W/(4.*eps_*eps_)
        W = W+W.T
        return W

def f_explore(fdyn,x0,time,eps_ = False):
        ''' compute the perturbation around a given trajectory  defined by the initial conditions x0.
        output has dimensions (nt x nx) x (nx x 2).
        (nt x nx) trajectories
        (nx x 2) positive and negative perturbation around each component.
        '''
        nx = x0.size
        nt = time.size
        
        dx = np.zeros((nt,nx,nx,2))
        if isinstance(fdyn,chemreac):#is fdyn a "cantera" object
                # for each species/comp of the state space
                for ix in range(nx):
                        #perturbation
                        base = np.zeros(nx)
                        if type(eps_) is float:base[ix] = eps_
                        else:base[ix]=eps_[ix]
                        #def new reactions and integrate them
                        reaction_p = chemreac(cti=fdyn.cti,V=fdyn.V,problem = fdyn.problem)
                        reaction_p.TP(fdyn.tpinit,reset=False)
                        reaction_p.concentrations(c=x0+base,reset=True)
                        reaction_p.integrate(time=time)
                        xp = np.array(reaction_p.h_z)[1:,:]
                        reaction_m = chemreac(cti=fdyn.cti,V=fdyn.V,problem = fdyn.problem)
                        reaction_m.TP(fdyn.tpinit,reset=False)
                        reaction_m.concentrations(c=x0-base,reset=True)
                        reaction_m.integrate(time=time)
                        xm = np.array(reaction_m.h_z)[1:,:]
                        # update the data
                        dx[:,:,ix,0] = xp
                        dx[:,:,ix,1] = xm
        else:# then it is a dynamical system
                for ix in range(nx):
                        base = np.zeros(nx)
                        base[ix] = eps_
                        xp = ode(fdyn,time,x0+base)
                        xm = ode(fdyn,time,x0-base)
                        dx[:,:,ix,0] = xp
                        dx[:,:,ix,1] = xm
        return dx
def f_explore_save(fdyn,x0,time,eps_ = False):
        ''' compute the perturbation around a given trajectory  defined by the initial conditions x0.
        output has dimensions (nt x nx) x (nx x 2).
        (nt x nx) trajectories
        (nx x 2) positive and negative perturbation around each component.
        '''
        nx = x0.size
        nt = time.size
        
        dx = np.zeros((nt,nx,nx,2))
        if isinstance(fdyn,chemreac):#is fdyn a "cantera" object
                # for each species/comp of the state space
                for i in range(nx):
                        #perturbation
                        base = np.zeros(nx)
                        base[i] = eps_
                        #def new reactions and integrate them
                        reaction_p = chemreac(cti=fdyn.cti,V=fdyn.V,problem = fdyn.problem,concentrations=x0+base)
                        reaction_p.integrate(time=time)
                        xp = np.array(reaction_p.h_z)[1:,:]
                        reaction_m = chemreac(cti=fdyn.cti,V=fdyn.V,problem = fdyn.problem,concentrations=x0-base)
                        reaction_m.integrate(time=time)
                        xm = np.array(reaction_m.h_z)[1:,:]
                        # update the data
                        dx[:,:,i,0] = xp
                        dx[:,:,i,1] = xm
        else:# then it is a dynamical system
                for i in range(nx):
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
        for i in range(nx):
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


def M(W,measure='trace'):
        ''' Measure of the quality of the Gramian '''
        if measure is 'cond':
                J = np.linalg.cond(W)
        if measure is 'det':
                try:
                        J=np.linalg.det( np.linalg.inv(W) )
                except:
                        J = 1000
        if measure is 'trace':
                try:
                        J=np.trace( np.linalg.inv(W) )
                        if J<0:J=np.abs(J)
                except:
                        J = 1000
        if measure is 'mixte':J=M(W,'trace')+M(W,'det')

        if J == np.inf:
                J = 1.e18
        
        return np.log(J)

def ode(fdyn,time,x0):
        '''Integrate system at hand with Runge Kutta 4'''
        dt = time[1]-time[0]
        nt = time.size
        nx = x0.size
        X = np.zeros((nt,nx))
        for it,t in enumerate(time):
                if it==0: X[it,:] = x0
                else: X[it,:] = dyn.rk4(fdyn,X[it-1,:],dt,t)
        return X
#########################################################
############ TO FIX/Not properly implented

def cost_neighboor(alpha,Wi, alpha_0):  
        W = Wi[0]
        for i in range(len(alpha)): W += alpha[i]*Wi[i+1]
        J = - M(W,'cond') + 0.*np.linalg.norm(alpha-alpha_0)
        # check tangent ?
        return J

def alpha_grid(bounds,n):
        alphas = np.ogrid[[slice(row[0], row[1], i_n*1j) for row, i_n in zip(bounds, n)]]
        #bounds = np.repeat([(0,1)], D, axis = 0)       
        return alphas

def grid2list(alphas):
        pass

def response_(alphas=None,dx=None,f_def_obs=None,time=None,eps_=None,rho=1.e-3):
        if alphas is None:
                print('please provide alphas')
        S=[]
        if dx is None:
                print('run f_explore first')
                return -1
        if time is None:
                print('time needed')
                return -1
        if eps_ is None:
                print('need precision')
                return -1
        if f_def_obs is None:
                print(' obserabl def as y=sum alpha_i x_i')
                f_def_obs = lambda x,alpha: np.dot(alpha,x)
        J = lambda alpha: cost_gramian(alpha, f_def_obs=f_def_obs, rho=rho, dx=dx, eps_=eps_, time=time)
        for alpha in alphas:
                S.append(J(alpha))
        return S

        dx0 = np.zeros((2,n,n))
