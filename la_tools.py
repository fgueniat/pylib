import numpy as np



def gradient2d(u,x,axis=0):
    if axis==1:
        du = np.c_[u[:,1]-u[:,0],(np.diff(u,axis=1)[:,1:]+np.diff(u,axis=1)[:,0:-1])/2,u[:,-1]-u[:,-2]]
        dx = np.c_[x[:,1]-x[:,0],(np.diff(x,axis=1)[:,1:]+np.diff(x,axis=1)[:,0:-1])/2,x[:,-1]-x[:,-2]]
    if axis==0:
        v=np.copy(u).T
        du = np.c_[v[:,1]-v[:,0],(np.diff(v,axis=1)[:,1:]+np.diff(v,axis=1)[:,0:-1])/2,v[:,-1]-v[:,-2]]
        del v
        du=du.T
        y=np.copy(x).T
        dx = np.c_[y[:,1]-y[:,0],(np.diff(y,axis=1)[:,1:]+np.diff(y,axis=1)[:,0:-1])/2,y[:,-1]-y[:,-2]]
        dx=dx.T
        del y
    dudx = np.divide(du,dx)
    return dudx


def curl(u,v,x=False,y=False):
	''' return the 2D vorticity field
	'''
	if x is False:x=np.ones(u.shape)
	if y is False:y=np.ones(u.shape)
	
        #		
        #	try:
        #		dx = np.c_[x[:,1]-x[:,0],(np.diff(x,axis=1)[:,1:]+np.diff(x,axis=1)[:,0:-1])/2,x[:,-1]-x[:,-2]] 
        #		dy = np.r_[[y[1,:]-y[0,:]],(np.diff(y,axis=0)[1:,:]+np.diff(y,axis=0)[0:-1,:])/2,[y[-1,:]-y[-2,:]]] 
        #	except:
        #		print 'try to use meshgrid'
        #		return -1
        #		
	#c = np.gradient(u,dy,axis = 0) - np.gradient(v,dx,axis=1)
	c = gradient2d(u,y,axis = 0) - gradient2d(v,x,axis=1)

	return c

def ssp2(f, g, t, x0, u, p):

	# Local Variables: xk, g, f, h, K, p, s, u, t, uk, hs, y, x, x0, STAGES, k, tk
	# Function calls: isscalar, ssp2, floor
	#%%% summary: ssp2 (Low-Storage Stability Preserving Runge-Kutta SSP32)
	#%$
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

def sparse(i=False,j=False,v=False,m=False,n=False):
	if i is False or j is False: 
		print('improper arguments')
		return -1
	else:
		i = np.atleast_1d(i)
		j = np.atleast_1d(j)
	if v is False:
		print('v is false')
		return np.zeros((i[0],j[0]))
	else:
		v = np.atleast_1d(v)
	if m is False: m = np.max(i)+1
	if n is False: n = np.max(j)+1
	S = np.zeros((m,n))
	for ind in range(len(i)):
		S[i[ind],j[ind]] = v[ind]
	return S

def mdot(ABC):
    ''' Multiply several matrix at the same time'''
    A = ABC[0]
    for i in range(1,len(ABC)):
    	A = np.dot(A,ABC[i])
    return A

def projector(v):
	return np.dot(np.linalg.pinv(np.dot(v.T,v)),v.T)

def RMSE(a,b):
	err=np.sqrt(  np.sum( (a-b)*(a-b) )/a.size  )
	return err

def gen_cholesky(H, eps_machine=1.0e-15, print_prefix=0, print_flag=0): 
#The Gill, Murray, and Wright modified Cholesky algorithm, from Practical Optimization, Academic Press, London, 1981
 
	ndim = len(H) 
	I = np.eye(ndim) 
	# Calculate gamma(A) and xi(A). 
	gamma = 0.0 
	xi = 0.0 
	for i in range(ndim): 
		gamma = np.max((np.abs(H[i, i]), gamma)) 
		for j in range(i+1, ndim): 
			xi = max(np.abs(H[i, j]), xi) 

# Identify delta and beta. 
	try:
		delta = eps_machine * np.max((gamma + xi, 1.0)) 
	except:
		print(gamma)
		print(xi)
		print(eps_machine)
		delta = eps_machine * np.max((gamma + xi, 1.0)) 

	if ndim == 1: 
		beta = np.sqrt( np.max((gamma, eps_machine)) )
	else: 
		beta = np.sqrt(  np.max((gamma, xi / np.sqrt(ndim**2 - 1.0), eps_machine))  ) 
 
# Initialise data structures. 
	a = 1.0 * H 
	r = 0.0 * H 
	P = 1.0 * I 

# Main loop. 
	for j in range(ndim): 
# Row and column swapping, find the index > j of the largest diagonal element. 
		q = j 
		for i in range(j+1, ndim): 
			if np.abs(a[i, i]) >= np.abs(a[q, q]): 
				q = i 

# swap row and column j and q (if j != q). 
		if q != j: 
# Temporary permutation matrix for swaping 2 rows or columns. 
			p = 1.0 * I 

# Modify the permutation matrix P by swaping columns. 
			row_P = 1.0*P[:, q] 
			P[:, q] = P[:, j] 
			P[:, j] = row_P 

# Modify the permutation matrix p by swaping rows (same as columns because p = pT). 
			row_p = 1.0*p[q] 
			p[q] = p[j] 
			p[j] = row_p 

# Permute a and r (p = pT). 
		a = np.dot(p, np.dot(a, p)) 
		r = np.dot(r, p) 
 
# Calculate dj. 
		theta_j = 0.0 
		if j < ndim-1: 
			for i in range(j+1, ndim): 
				theta_j = np.max(  (theta_j, np.abs(a[j, i]))  ) 
		dj = np.max(  (np.abs(a[j, j]), (theta_j/beta)**2, delta)  ) 

# Calculate row j of r and update a. 
		r[j, j] = np.sqrt(dj)     # Damned sqrt introduces roundoff error. 
		for i in range(j+1, ndim): 
			r[j, i] = a[j, i] / r[j, j] 
			for k in range(j+1, i+1): 
				a[i, k] = a[k, i] = a[k, i] - r[j, i] * r[j, k]     # Keep matrix a symmetric. 
 
# Finally, the Cholesky decomposition of H. 
	L = np.dot(P, np.transpose(r))
	return L


def Smooth_random(n,x,m,var):
	noise = np.zeros(n)
	if np.abs(var)>1.e-12:
		for i in range(1,m):
			noise += np.random.normal(0,np.sqrt(var)) * np.sin(i * np.pi * x)
			if np.mod(i,2)==1:
				noise += np.random.normal(0,np.sqrt(var)) * np.cos(i * np.pi * x / 2.0)
	return noise

def Smooth_cov(n,x,m,var,isL=False):
	L = np.zeros((m,n))
	for j in range(0,m):
		for k in range(0,n):
			if np.mod(k,2) ==1:
				L[j,k] = np.sqrt(var) * (  np.sin(np.pi * j * x[k]) + np.cos(np.pi * j * x[k])  )
			else:
				L[j,k] = np.sqrt(var) * np.sin(np.pi * j * x[k])
	if isL is True:
		return L
	else:
		return np.dot(L.T,L)

def MG2D(x,lissage=1,filtre='lin'):
    ''' mean filter in 2d
    '''
    def mg(x,lissage,filtre):
        xl = np.empty(x.shape)
        for ix,xx in enumerate(x):
            xl[ix,:] = MG(xx,lissage=lissage,filtre=filtre)
        return xl
    xl = mg(x,lissage=lissage,filtre=filtre)
    xl = mg(xl.T,lissage=lissage,filtre=filtre).T
    return xl



def smooth(x,lissage=1,filtre = "lin",recursion = 1):
    ''' binding to MG'''
    return MG(x,lissage=lissage,filtre=filtre,recursion=recursion)

def MG(x,lissage=1,filtre = "lin",recursion = 1):
    '''
    Mean filtering
    '''
    if recursion == 0:
        return x
    #
    X = np.copy(x)
    y = np.zeros(x.size)
    lin = lambda n: 1.0*np.ones(2*n+1) / (2*n+1)
    wma = lambda n: 1.0* np.append(np.arange(1,n+1),np.arange(n+1,0,-1)) / (((n+1)*(n+2))/2 + ((n+1)*(n))/2)
    # not implemented, wma is used
    expo = lambda n: wma(n)
    while recursion>0:
        recursion = recursion-1        
        for i in range(x.size):
            if i<lissage:ni = i
            elif x.size-i -1 < lissage: ni = x.size-i-1
            else: ni = lissage
            if filtre == "lin": A = lin(ni)
            elif filtre == "wma": A = wma(ni)
            elif filtre == "expo": A = expo(ni)
            else: A = lin(ni)
            y[i] = np.dot(A,X[i-ni:i+ni+1])
        X = np.copy(y)

    return y

def MGfast(x,lissage=1):
	
	return x

def MedG(x,lissage=1):
	y = np.zeros(x.size)

	for i in range(x.size):
		if i<lissage:
			ni = i
		elif x.size-i -1 < lissage:
			ni = x.size-i-1
		else:
			ni = lissage
		y[i] = np.median(x[i-ni:i+ni+1])
	return y

def MMG(x,lissage=1,filtre = "lin"):
	y = MedG(x,lissage)
	y = MG(y,lissage,filtre)
	return y
def FG(x,lissage=0.9):
    ''' 
    Fourier Filtering
    '''
    from scipy import fft,ifft
    yfft=fft(x)
    for i in range(x.size):
            if i<x.size/2.0:
                    n = 2*i/x.size
            else:
                    n = -2*(x.size/2.0 - i)/x.size
     #       print n
            if n>lissage:
                    yfft[i] = 0.0
    #print(yfft)
    y=ifft(yfft)
    #print(y-x)
    return np.real(y)

def clear(x):
    y = np.array(x)
    y =y.flatten()
    y = y[~np.isnan(y)]
    y = y[~np.isinf(y)]
    return y
def mean(x,return_x = False):return mean_outliers(x,return_x)
def mean_outliers(x,return_x = False):
    '''from Leys, Christophe, et al. "Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median." Journal of Experimental Social Psychology 49.4 (2013): 764-766.'''
    b = 1.4826
    try:
            y = clear(x)
    except:
            return 0.
    if y.size == 0:
            return 0
    mx = np.median(y)
    mad = b*np.median(np.abs(y-mx))
    if mad < 1.e-16:return np.mean(y)
    x2 = y[((y<mx+2.0*mad)&(y>mx-2.0*mad))]
    if return_x is True:
        return x2
    return np.nanmean(x2)

def std(x):return std_outliers(x)
def std_outliers(x):
    '''from Leys, Christophe, et al. "Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median." Journal of Experimental Social Psychology 49.4 (2013): 764-766.'''
    b = 1.4826
    try:
        y = clear(x)
    except:
        return 0.
    if y.size == 0:
        return 0
    mx = np.median(y)
    mad = b*np.median(np.abs(y-mx))
    if mad < 1.e-16:mad = b*mean_outliers(np.abs(y-mx))
    if mad < 1.e-16:return 0.
    x2 = y[((y<mx+2.0*mad)&(y>mx-2.0*mad))]
    stdx = np.std(x2)
    if stdx < 1.e-16 and np.linalg.norm(x2-mx)/np.linalg.norm(y)<1.e-16:
        return np.std(y)
    else:
        return stdx

def find_non_outliers(x):
    '''from Leys, Christophe, et al. "Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median." Journal of Experimental Social Psychology 49.4 (2013): 764-766.'''
    b = 1.4826
    try:
        y = clear(x)
    except:
        return 0.
    if y.size == 0:
        return 0
    mx = np.median(y)
    mad = b*np.median(np.abs(y-mx))
    if mad < 1.e-16:mad = b*mean_outliers(np.abs(y-mx))
    if mad < 1.e-16:return 0.
    inds = np.where((y<mx+2.0*mad)&(y>mx-2.0*mad))
    return inds
