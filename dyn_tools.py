import numpy as np
import scipy as sp
from la_tools import projector

#mport weave

def delta_(t,h):
#% Delta Impulse Input
        d=np.copy(np.atleast_1d(t))
        d[d>h] = 0.
        d[d<0.] = 0.
        return d
 
def DxDt_r(x ,t ,p=False ):
        ''''Roessler '''
        if p is False: a,b,c = 0.398,2.,4.0
        else: a,b,c = p[0],p[1],p[2]
        dx = np.zeros(3)
        dx[0] = - x[1] - x[2]
        dx[1] = x[0] + a*x[1]
        dx[2] = b + x[2] *(x[0] - c)
        return dx


def DxDt_mass(x ,t ,p=False ):
        ''''double dumped mass system '''
        w2 = 1.0
        dx = np.zeros(4)
        dx[0] =  x[2] 
        dx[1] = x[3]
        dx[2] = -2*w2*x[0] + w2*x[1]
        dx[3] = w2*x[0] - 2*w2*x[1]
        return dx


def f_ds(x ,t ,p=False ):
        a,b = -8./7,-5./7
        return b*x + 0.5 * (a - b) * (np.abs(x+1) - np.abs(x-1) )

def DxDt_ds(x ,t ,p=False ):
        ''''Double scroll '''
        if p is False:alpha, beta, gamma = 9.,100./7,0.
        else: alpha, beta, gamma = p[0],p[1],p[2]
        dx = np.zeros(3)
        dx[0] = alpha * ( x[1] - x[0] - f_ds(x[0]) )
        dx[1] = x[0] - x[1] + x[2]
        dx[2] = -beta * x[1] -gamma * x[2]
        return dx

def DxDt_l(x ,t ,p =False):
        ''''Lorenz '''
        if p is False: sigma, rho, beta = 10., 28., 8./3.
        else: sigma, rho, beta = p[0],p[1],p[2]
        dx = np.zeros(3)
        dx[0] = sigma* (x[1] - x[0])
        dx[1] = x[0]*( rho-x[2])  - x[1]
        dx[2] = x[0] * x[1] - beta* x[2]
        return dx

def DxDt_cord(x ,t ,p =False):
        ''''Cord attractor'''
        if p is False: a,b,F,G = 0.25,4.0,8.0,1.0       
        else: a,b,F,G = p[0],p[1],p[2], p[3]    
        dx = np.zeros(3)
        dx[0] = -x[1] -x[2] - a*x[0] + a*F
        dx[1] = x[0]*x[1]  - b*x[0]*x[2] - x[1] + G
        dx[2] = b*x[0]*x[1] + x[0]*x[2] - x[2]
        return dx

def DxDt_l84(x ,t ,p =False):
        ''''Lorenz 84'''
        if p is False: a,b,F,G = 0.25,4.0,8.0,1.0       
        else: a,b,F,G = p[0],p[1],p[2], p[3]    
        dx = np.zeros(3)
        dx[0] = -x[1]*x[1] -x[2]*x[2] - a*x[0] + a*F
        dx[1] = x[0]*x[1]  - b*x[0]*x[2] - x[1] + G
        dx[2] = b*x[0]*x[1] + x[0]*x[2] - x[2]
        return dx


def rk4(f,x,dt,t=False,p=False,u=False,scipy=False):
    ''' rk4 integrator
def rk4(f,x,dt,t=False,p=False,u=False)'''
    if scipy is True:
        if t is False: t=0.
        if u is False: ut= lambda t: 0.
        if callable(u) is True: ut = lambda t: u(t)
        else: ut= lambda t: u
        if p is False: 
            fdyn = lambda x,t: f(x,t)
        else: fdyn = lambda x,t,p: f(x,t,p)
        return sp.integrate.odeint(fdyn,x,[t,t+dt])[-1,:]
    else:
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

def rk4_c(f,x,dt,t=False,p=False,u=False):
    ''' rk4 integrator using weave
def rk4(f,x,dt,t=False,p=False,u=False)'''
        # not implemented yet
    pass

def integrate(f,x0,dt, t0=0, tfin = 10., p=False,u=False):
        y=[]
        n = int((tfin-t0)/dt)
        time = t0+np.arange(n)*dt
        for it, t in enumerate(time):
                if it==0:
                    if isinstance(x0,np.ndarray): y.append(x0.tolist())
                    else: y.append(x0)
                else:y.append(rk4(f,y[-1],dt,t,p,u))
        return y,time.tolist()

def embedding(x,d_emb = 3,type_emb=None,dt = False,n_emb = 10):
    '''
    type_emb: 
    * None = time difference, with n_emb
    * 'svd': time differnece, then rescalde with svd
    * 'deriv': numerical time deriv using dt. Not implemented yet
    '''
    if type_emb is not 'deriv':
        X = np.array(
            [x[0+i*n_emb:-(d_emb-i)*n_emb] for i in range(d_emb)]
            ).T
    else:
        print('not implemented yet')
        return -1
    if type_emb is 'svd':
        u,s,v = np.linalg.svd(X,0)
        P = projector(v)
        X = np.dot(X,P)
    return X
        

def visibility_graph(x,time=False,n=20,mem_save = False,forced=False,verbose = False,tresh = 100,fast = True):
    '''
    visibility graph for discriminating chaos and stochastic
    PHYSICAL REVIEW E 80, 046103 2009
    Horizontal visibility graphs: Exact results for random time series
    B. Luque,L. Lacasa, F. Ballesteros and J. Luque
    '''
    if fast is True: fast = False #deprecated
    if time is False:
        time = np.arange(len(x))
    if len(time) != len(x):
        print('lengths of time and time series have to match')
        return -1
    
    nx,nt = len(x),len(time)
    
    #rescale
    X = x-np.min(x)
    X = X/np.max(X)
    #
    if nx>1500:
        if forced is not True:
            prin('mem_save turned on')
            mem_save = True

    if mem_save is False:
        graph = np.zeros((nx,nx))
        print_prog = 0
        last_print = 0
        nmax = nt
        for ix in range(nx-2):
            for iy in range(ix+2,nx):
                if verbose is True:
                    print_prog +=1
                    if int(100.*print_prog/nmax) > last_print:
                        last_print +=1
                        print('graph: ' + str(last_print -1) + '% computed')
                if ix==iy:
                    pass
                else:
                    ya = X[ix]
                    yb = X[iy]
                    ta = time[ix]
                    tb = time[iy]
                    test = True
                    for ic in range(ix+1,iy):
                        yc = X[ic]
                        tc = time[ic]
                        if yc > ya + (yb-ya)*(tc-ta)/(tb-ta):
                            test = False
                            break
                    if test is False:
                        graph[ix,iy] = 0
                        graph[iy,ix] = 0
                    else:
                        graph[ix,iy] = 1
                        graph[iy,ix] = 1
    else:        
        print_prog = 0
        last_print = 0
        nmax = nt
        graph = np.zeros(nx)
        for ix in range(nx-1):            
            if verbose is True:
                    print_prog +=1
                    if int(100.*print_prog/nmax) > last_print:
                        last_print +=1
                        print('graph: ' + str(last_print -1) + '% computed')
            
            if fast is True:
                #if  (yc > ya + (yb-ya)*(tc-ta)/(tb-ta)) 
                code = '''
                    double ya,yb,yc;
                    double ta,tb,tc;
                    int test;
                    for (int iy=ix+1;iy<nx;++iy){
                        ya = X(ix);
                        yb = X(iy);
                        ta = time(ix);
                        tb = time(iy);
                        test = 1;
                        for (int ic = ix+1;ic<iy;++ic){
                            yc = X(ic);
                            tc = time(ic);
                            if  (yc > ya) {
                                test = 0;
                                break;
                                }
                            if  (yc > yb) {
                                test = 0;
                                break;
                                }
                            }
                        if (test == 1){
                            graph(ix) +=1;
                            graph(iy) +=1;
                            }
                        }

                    '''
                #weave.inline(code,['graph','X','time','nx','ix'],type_converters=weave.converters.blitz)
            else:
                for iy in range(ix+2,nx):
                    if ix-iy<tresh:
                        if ix==iy:
                            pass
                        else:
                            ya = X[ix]
                            yb = X[iy]
                            ta = time[ix]
                            tb = time[iy]
                            test = True
                            for ic in range(ix+1,iy):
                                yc = X[ic]
                                tc = time[ic]
                                if yc > ya + (yb-ya)*(tc-ta)/(tb-ta):
                                    test = False
                                    break
                            if test is False:
                                pass
                            else:
                                graph[ix] += 1
                                graph[iy] += 1
    return graph



def graph_distrib(graph,plot=False,kmax = None):
    if kmax is None:
        k = np.arange(np.max(graph)+1)
    else:
        k = np.arange(kmax+1)
    pk = np.zeros(len(k))
    for g in graph:
        if kmax is None:
            pk[g]+=1
        else:
            if g<kmax+1:
                pk[g]+=1

    k=k[1:]
    pk = pk[1:]
    
    pk = 1.*pk/np.sum(pk)
    pk[pk==0.] = np.min(pk[np.nonzero(pk)])/10.
    
    A = np.ones((len(k)-2,2))
    A[:,1]=k[1:-1]
    logpk = np.log(pk[1:-1])
    mloga,lambda_ = -np.linalg.lstsq(A,logpk)[0][:]
    a = np.exp(-mloga)
    print('lambda approx', lambda_)
    if plot is True:
        import plot_tools
        plot_tools.multiplot2((k,k),(pk,a*np.exp(-lambda_*k)),plot_tools.Paradraw(marks = ['-','--'], markers=['o',''],y_scale = 'log',x_label = 'k',y_label = 'P(k)',ax_x_format = '%.0f',legend = ['','$\lambda = ' + '{0:.2f}'.format(lambda_) + '$']))
    return k,pk,lambda_

def Poincarre(s,r=0.):
    time_series = s-r
    dts = np.diff(time_series)
    pass


