import numpy as np
from scipy.optimize import minimize as solve
from scipy.optimize import basinhopping as recuit
from scipy.linalg import sqrtm
import dyn_tools as dyn
import plot_tools as pt
from la_tools import mean_outliers as smart_mean
from cantera_tools import chemreac
#import weave

import time as temps

def params(argv=None,inputs = None):
    '''
    argv: parameters can be sent directly through the console
        argv #1: int. it corresponds to the constrains associated to the [int] species
        argv #2: int. Number of linear constrains
        argv #3: int. Number of quadratic constrains
        argv #4: int. Number of cubic constrains
        argv #5: tuple. initial concentrations
        argv #6: tuple. initial TP
        argv #7: tuple. Second set of initial concentrations, for comparison purpose
        argv #8: tuple. Second set of initial TP, for comparison purpose


    Parameters for the ID, as a dictionary
            'input_mechanism' cti file e.g., 'boivin.cti'
            'input_problem' type of problem, e.g. ('0D','pressure'),
            'input_V' initial volume
            'input_nt' number of step size, e.g. 150
            'input_tmin' start time, e.g. 1.e-7,
            'input_tmax' end time, e.g. 1.e-2
            'input_offset_ts' offset for focus on given time scale. Has to be between 0 and 1, e.g. 0.2
            'input_offtype' set to 'slow' or 'fast', to focus on different time scales
            'input_eps_' parameter for the gramian. 5.e-2
            'input_rho' parameter of the gramian 1.e-9,
            'input_measure' type of metric for the gramian, e.g. 'trace'
            'input_spec' input species. its constrains associated to the slow manifold will be identified, e.g.'O2',
            'input_method' method used for the minimization, e.g. 'BFGS',
            'input_index_sorted' None,
            'input_n' int. Number of linear constrains
            'input_n_square' int. Number of quadratic constrains
            'input_n_cubic' int. Number of cubic constrains
            'input_tp' tuple. initial TP
            'input_ci' tuple. initial concentrations


            'input_ci_random' Second set of initial concentrations, for comparison purpose. Default is None
            'input_tp_random' Second set of initial TP, for comparison purpose. Default is None
            'input_minimization' preexisting results for avoiding re-computations. Default is None.
    '''
    if inputs is None:
            inputs = {
            'input_mechanism':'boivin.cti',
            'input_problem':('0D','pressure'),
            'input_V':1.e3 * (1.e-2)**3,
            'input_nt':150,
            'input_tmin':1.e-7,
            'input_tmax':1.e-2,
            'input_offset_ts':0.2,
            'input_offtype':'slow',
            'input_eps_':5.e-2,
            'input_rho':1.e-9,
            'input_measure':'trace',
            'input_spec':'O2',
            'input_method':'BFGS',
            'input_index_sorted':None,
            'input_n':4,
            'input_n_square':0,
            'input_n_cubic':0,
            'input_tp':(1200,101325),
            'input_ci':np.array([1.e-03,1.0e-09,1.e-03,1.e-09,1.e-9,1.e-09,1.e-9,1.e-9]),
            'input_ci_random':None,
            'input_tp_random':None,
            'input_minimization':None,
            }
    try:
        inputs['input_spec'] = argv[1]
        inputs['input_n'] = int(argv[2])
        inputs['input_n_square'] = int(argv[3])
        inputs['input_n_cubic'] = int(argv[4])
        inputs['input_ci'] = dt.split_list(argv[5])
        inputs['input_tp'] = dt.split_tuple(argv[6])
        inputs['input_ci_random'] = dt.split_list(argv[7])
        inputs['input_tp_random'] = dt.split_tuple(argv[8])
    except:
        pass    
    inputs['spec'] = inputs['input_spec']
    return inputs





measures_possible = ['cond', 'det','trace','max_svd']

def slow_manifold(alpha_0, dx=None, f_def_obs=None,eps_ = None,time = None,rho = 1.e-3, measure = None, verbose = False,offset = 0.2,offtype='slow',method = 'Powell',isrecuit=True, niter=10,maxiter=[10,1],n_recuit=100,r=None):
    '''
    Identify geometric constrains based on observability analysis
        alpha_0 array. inital conditions for the minimization. Size is the number of constrains.
        dx array. Precomputed arrays for the gramian. Output from f_explore.
        f_def_obs function that define the observable. obs = f_def_obs(x) = sum a_i x_i + sum b_i x_i**2 + sum c_i x_i**3
        eps_ real. parameter for the gramian
        time array. times corresponding to the snapshots
        rho real. parameter for the minimization
        measure str. metric for the gramian, such as 'trace' 
        verbose bool.
        offset real. offset for focus on given time scale. Has to be between 0 and 1, e.g. 0.2
        offtype str. time scale of interst, e.g. 'slow'
        method str. minimzation method, e.g. 'Powell' or 'BFGS'
        isrecuit bool. annealing
        niter int. number of anneling
        maxiter arrat. number of iteration and function evaluation.
        n_recuit int. number of annealing
    
    '''
    #verbose = True
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
    argums = (f_def_obs,rho,dx,eps_,time,False,False,False,measure,offset,offtype,False,r)
    if type(maxiter) is not list: 
        if type(maxiter) is int or type(maxiter) is float:
            maxiter=[maxiter]
        else:
            print('maxiter ill defined')
            return -1
    if isrecuit is True:#Annealing
        opt = {'maxiter':maxiter[1],'disp':True,'maxfev':10}
        opt = {'maxiter':0,'disp':True,'maxfev':len(alpha_0)}
        print(opt)
        min_kwargs = {'args':argums,'method':'Powell','options':opt}        

        np.random.seed(111111)
        alpha_0 = recuit(cost_gramian,alpha_0,niter=n_recuit,minimizer_kwargs = min_kwargs,T=20.,disp=True)
        np.random.seed()
        print(alpha_0)
        alpha_0 = alpha_0.x
    opt = {'maxiter':maxiter[0],'disp':True}
    #opt = {'maxiter':5,'disp':True}
    #if no annealing, we finish with a regular minimization
    try:
        a = solve(cost_gramian,alpha_0,method = method, tol = 1.e-5,args = argums,options=opt)
    except:
        a = solve(cost_gramian,alpha_0,method = 'BFGS', tol = 1.e-5,args = argums,options=opt)

    return a

def cost_gramian(alpha, f_def_obs=False, rho=0., dx=False, eps_=False, time=False, nt = False, ny=False, nx = False, measure=None, offset = 0.2, offtype ='slow', verbose = False,r=None):
    '''
    Compute the gramian and its  measure
        alpha array. inital conditions for the minimization. Size is the number of constrains.
        f_def_obs function that define the observable. obs = f_def_obs(x) = sum a_i x_i + sum b_i x_i**2 + sum c_i x_i**3
        rho real. parameter for the minimization
        dx array. Precomputed arrays for the gramian. Output from f_explore.
        eps_ real. parameter for the gramian
        time array. times corresponding to the snapshots
        nt int. dimension of time
        ny int. dimension of observable
        nx int. dimension of state vector
        measure str. metric for the gramian, such as 'trace' 
        offset real. offset for focus on given time scale. Has to be between 0 and 1, e.g. 0.2
        offtype str. time scale of interst, e.g. 'slow'
        verbose bool.
    '''
#    print np.random.randint(0,10)
#    print int(temps.time())
    # parameters
    if f_def_obs is False:#def obs
        fobs = lambda xx: xx[1] - np.sum([xx[0]**a for a in alpha])
    else:fobs = lambda x: f_def_obs(x,alpha)
    if nt is False:nt = np.size(time)
    if nx is False:nx = np.size(dx[0,:,0,0])
    if ny is False:ny = np.size(fobs(dx[0,:,0,0]))
    #
    #computation of the gramian
    w = EG(fobs = fobs, dx=dx, eps_ = eps_, time = time, nt = nt, nx = nx, ny = ny, offset = offset,r=r,offtype = offtype)
    if measure is None or measure not in measures_possible:
        measure = 'trace'
    norm_type = 1    
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


def EG(fobs,dx,eps_,time,nt = False,nx = False,ny = False,offset = .2,r=None,offtype = 'slow',h_r=None):
    ''' 
    Compute the Empirical gramian
        fobs function that define the observable. obs = f_def_obs(x) = sum a_i x_i + sum b_i x_i**2 + sum c_i x_i**3
        dx array. Precomputed arrays for the gramian. Output from f_explore.
        eps_ real. parameter for the gramian
        time array. times corresponding to the snapshots
        nt int. dimension of time
        ny int. dimension of observable
        nx int. dimension of state vector
        offset real. offset for focus on given time scale. Has to be between 0 and 1, e.g. 0.2
        offtype str. time scale of interst, e.g. 'slow'
        verbose bool.  
    '''
    fast = False
    # parameters
    if nt is False:nt = np.size(time)
    if nx is False:nx = dx[0,:,0,0].size
    if ny is False:ny = fobs(dx[0,:,0,0]).size
    if isinstance(eps_, (list, tuple, np.ndarray)) is True:
        if len(eps_)<nx:
            epsil = eps_[0]* np.ones(nx)
        else:
            epsil = eps_* np.ones(nx)
    else:
        epsil = eps_* np.ones(nx)

    offset = int(offset*nt)+1
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
    if fast is False:
        for it in range(nstart,nend):
            dt = time[it]-time[it-1]
            y=np.zeros((nx,2*ny))
            if r is None:
                if h_r is not None:
                    r0=h_r[it]
                else:
                    r0 = np.zeros(nx)
            else: r0 = r.h_z[it]
            for ix in range(nx):# observation
                epsi=2*epsil[ix]
                y[ix,0:ny] = (fobs(dx[it,:,ix,0] - r0).flatten())/epsi
                y[ix,ny:] = (fobs(dx[it,:,ix,1] - r0).flatten())/epsi
            # contribution of each observation to the gramian    
            W = W+dt*np.dot(y,y.T) # equiv to W = W + y*y.T
    else:
        f = lambda x,e: fobs(x).flatten()/e
        dx_temp = np.zeros(dx[0,:,0,0].shape)
        for it in range(nstart,nend):
            dt = time[it]-time[it-1]
            y=np.zeros((nx,2*ny))
            if r is None:r0=np.zeros(nx)
            else: r0 = r.h_z[it]
            code = ''' 
            py::tuple arg(2);
            PyObject *W_temp;
            PyObject *y_temp;
            double *y_iy;
            double *y_il;
            double *y_jl;
            double *w_ij;
            
            W_temp = (PyObject*) PyArray_FromAny( W, PyArray_DescrFromType(NPY_FLOAT64),0,0,0,NULL);

            for (int ix=0;ix<nx;++ix){
                
                // fill dx_temp
                for (int i_dx = 0;i_dx < PyArray_DIM(dx_temp,0);++i_dx){
                    dx_temp(i_dx) = dx(it,i_dx,ix,0) - r0(i_dx); 
                    }
                arg[0] = dx_temp;
                arg[1] = 2 * epsil(ix);
                y_temp = (PyObject*) PyArray_FromAny( f.call(arg), PyArray_DescrFromType(NPY_FLOAT64),0,0,0,NULL);
                // fill y
                for (int i_y = 0;i_y < ny;++i_y){
                    y_iy = (double*) PyArray_GETPTR1(y_temp,i_y);
                    y(ix,i_y) = *y_iy;
                    }

                // fill dx_temp
                for (int i_dx = 0;i_dx < PyArray_DIM(dx_temp,0);++i_dx){
                    dx_temp(i_dx) = dx(it,i_dx,ix,1) - r0(i_dx); 
                    }
                arg[0] = dx_temp;
                arg[1] = 2 * epsil(ix);
                y_temp = (PyObject*) PyArray_FromAny( f.call(arg), PyArray_DescrFromType(NPY_FLOAT64),0,0,0,NULL);
                // fill y
                for (int i_y = 0;i_y < ny;++i_y){
                    y_iy = (double*) PyArray_GETPTR1(y_temp,i_y);
                    y(ix,ny+i_y) = *y_iy;
                    }

            // y is filled
            // now W +=  dt*y*y.T
            
            for (int i_i =0; i_i < nx ; ++i_i){
                for (int i_j =0; i_j < nx ; ++i_j){
                    w_ij = (double*) PyArray_GETPTR2(W_temp,i_i,i_j);
                    for (int i_l =0; i_l < 2 * ny ; ++i_l){
                        y_il = (double*) PyArray_GETPTR2(y,i_i,i_l);
                        y_jl = (double*) PyArray_GETPTR2(y,i_j,i_l);

                        *w_ij += dt * *y_il * *y_jl;
                        }
                    }
                }


            }
            return_val = W_temp;
          
            '''
            #W = weave.inline(code,['W','f','nx','ny','dx','dx_temp','r0','epsil','y'])
    #re normalization
    #W = W/(4.*eps_*eps_)
    W = W+W.T
    return W

def f_explore(fdyn,x0,time,eps_ = False):
    ''' 
    Compute the perturbation around a given trajectory defined by the initial conditions x0.
    Output has dimensions (nt x nx) x (nx x 2).
    (nt x nx) trajectories
    (nx x 2) positive and negative perturbation around each component.
        fdyn: function to integrate
        x0 array. intial conditions
        time array. times to get the snapshots
        eps_ array or real. perturbations around x0
        
    '''
    nx = x0.size
    nt = time.size
    
    dx = np.zeros((nt,nx,nx,2))
    eps_effectif = np.zeros(nx)
    if isinstance(fdyn,chemreac):#is fdyn a "cantera" object
        # for each species/comp of the state space
        for ix in range(nx):
            #perturbation
            base = np.zeros(nx)


            if isinstance(eps_, (list, tuple, np.ndarray)) is True:
                if len(eps_)<nx:
                    base[ix] = eps_[0]* x0[ix]
                else:
                    base[ix] = eps_[ix]
            else:
                base[ix] = eps_*x0[ix]

            #if type(eps_) is float:base[ix] = eps_*x0[ix]
            #else:base[ix]=eps_[ix]
            #def new reactions and integrate them
            reaction_p = chemreac(cti=fdyn.cti,V=fdyn.V,problem = fdyn.problem)
            reaction_p.TP(fdyn.tpinit,reset=False)
            reaction_p.Y(y=x0+base,reset=True)
            reaction_p.integrate(time=time)
            xp = np.array(reaction_p.h_z)[1:,:]
            reaction_m = chemreac(cti=fdyn.cti,V=fdyn.V,problem = fdyn.problem)
            reaction_m.TP(fdyn.tpinit,reset=False)
            reaction_m.Y(y=x0-base,reset=True)
            reaction_m.integrate(time=time)
            xm = np.array(reaction_m.h_z)[1:,:]
            # update the data
            dx[:,:,ix,0] = xp
            dx[:,:,ix,1] = xm
            #eps_effectif[ix] = (reaction_p.h_z[0][ix] - reaction_m.h_z[0][ix])/reaction_p.h_z[0][ix] 
            eps_effectif[ix] = reaction_p.h_z[0][ix] - reaction_m.h_z[0][ix] 
    else:# then it is a dynamical system
        for ix in range(nx):
            base = np.zeros(nx)
            if isinstance(eps_, (list, tuple, np.ndarray)) is True:
                if len(eps_)<nx:
                    base[ix] = eps_[0]* x0[ix]
                else:
                    base[ix] = eps_[ix]
            else:
                base[ix] = eps_ * x0[ix]
            #base[ix] = eps_
            xp = ode(fdyn,time,x0+base)
            xm = ode(fdyn,time,x0-base)
            dx[:,:,ix,0] = xp
            dx[:,:,ix,1] = xm
            eps_effectif[ix] = base[ix]
    return dx,eps_effectif
def f_explore_save(fdyn,x0,time,eps_ = False):
    ''' 
    Compute the perturbation around a given trajectory  defined by the initial conditions x0.
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
            reaction_p = chemreac(cti=fdyn.cti,V=fdyn.V,problem = fdyn.problem,Y=x0+base)
            reaction_p.integrate(time=time)
            xp = np.array(reaction_p.h_z)[1:,:]
            reaction_m = chemreac(cti=fdyn.cti,V=fdyn.V,problem = fdyn.problem,Y=x0-base)
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
    '''
    Simply plot the results of f_explore
    '''
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
    pt.multiplot2(xx,yy,pd)
    return xx,yy


def M(W,measure='max_svd'):
    ''' 
    Measure of the quality of the Gramian W
        W matrix. Gramian.
        measure str. type of metric 
    '''
    if measure is 'cond':
        J = np.linalg.cond(W)
    if measure is 'det':
        try:
            J=np.linalg.det( np.linalg.pinv(W,rcond=1.e-2) )
        except:
            J = 1000
    if measure is 'trace':
        try:
            J=np.trace( np.linalg.pinv(W,rcond=1.e-2) )
            if J<0:J=np.abs(J)
        except:
            J = 1000
    if measure is 'mixte':J=M(W,'trace')+M(W,'det')
    if measure is 'max_svd':
        u,s,v = np.linalg.svd(W)
        try:
            J=np.abs(np.log10(s[0])/np.log10(s[2]))
        except:
            J=np.abs(np.log10(s[0])/np.log10(s[1]))

    if J == np.inf:
        J = 1.e18
#    return J
    return J

def ode(fdyn,time,x0):
    '''
    Integrate system at hand with Runge Kutta 4
    '''
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




    
def print_obs(gram,specs,nullify = False, isplot=False,perm=True,fontsize = 14):
    ''' Plot a rescale of the observability properties'''
    if nullify is False: nullify = []
    W = np.sum(gram,axis=1)
    #
    W = np.log(W)
    W = W - np.median(W)
    for i_n in nullify:W[i_n] = 0
    W=W/np.sum(W)
    W+=-np.min(W)
    W=W/np.sum(W)
    for i_n in nullify:W[i_n] = 0
    s = np.std(W)
    m = smart_mean(W)
    W2 = np.zeros(W.shape)
    for iw,w in enumerate(W):
            if w < m+2*s:W2[iw]=w
    W2 = W2/np.sum(W2)
    print(W2)
    wm = np.max(W2)
    wM = np.max(W)
    for iw,w in enumerate(W):
            if w > m+2*s:
                    W2[iw]=wm*2*w/wM
    wm = np.min(W2[W2>0])
    for iw,w in enumerate(W2):
        if w <wm:
            if iw not in nullify: W2[iw] = 0.5*wm
    W = W2/np.sum(W2)

    if perm is False:
        perm = np.arange(len(specs))
    else:
        try:
            if len(perm) != len(specs):
                perm=np.argsort(W)[::-1]
        except:
                perm=np.argsort(W)[::-1]
    #
    if isplot is True:pt.plot1(np.array([W[pi] for pi in perm]),pt.Paradraw(x_tick_label = [specs[pi] for pi in perm],x_label = 'species',y_label = 'W',stem=True,ticksize = fontsize,thickness=['3']))
    return W

