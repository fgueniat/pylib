import numpy as np

def mdot(ABC):
	A = ABC[0]
	for i in xrange(1,len(ABC)):
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
	for i in xrange(ndim): 
		gamma = np.max((np.abs(H[i, i]), gamma)) 
		for j in xrange(i+1, ndim): 
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
	for j in xrange(ndim): 
# Row and column swapping, find the index > j of the largest diagonal element. 
		q = j 
		for i in xrange(j+1, ndim): 
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
			for i in xrange(j+1, ndim): 
				theta_j = np.max(  (theta_j, np.abs(a[j, i]))  ) 
		dj = np.max(  (np.abs(a[j, j]), (theta_j/beta)**2, delta)  ) 

# Calculate row j of r and update a. 
		r[j, j] = np.sqrt(dj)     # Damned sqrt introduces roundoff error. 
		for i in xrange(j+1, ndim): 
			r[j, i] = a[j, i] / r[j, j] 
			for k in xrange(j+1, i+1): 
				a[i, k] = a[k, i] = a[k, i] - r[j, i] * r[j, k]     # Keep matrix a symmetric. 
 
# Finally, the Cholesky decomposition of H. 
	L = np.dot(P, np.transpose(r))
	return L


def Smooth_random(n,x,m,var):
	noise = np.zeros(n)
	if np.abs(var)>1.e-12:
		for i in xrange(1,m):
			noise += np.random.normal(0,np.sqrt(var)) * np.sin(i * np.pi * x)
			if np.mod(i,2)==1:
				noise += np.random.normal(0,np.sqrt(var)) * np.cos(i * np.pi * x / 2.0)
	return noise

def Smooth_cov(n,x,m,var,isL=False):
	L = np.zeros((m,n))
	for j in xrange(0,m):
		for k in xrange(0,n):
			if np.mod(k,2) ==1:
				L[j,k] = np.sqrt(var) * (  np.sin(np.pi * j * x[k]) + np.cos(np.pi * j * x[k])  )
			else:
				L[j,k] = np.sqrt(var) * np.sin(np.pi * j * x[k])
	if isL is True:
		return L
	else:
		return np.dot(L.T,L)

def MG(x,lissage=1,filtre = "lin"):
	y = np.zeros(x.size)

	def lin(n):
		return 1.0*np.ones(2*n+1) / (2*n+1)
	def wma(n):
		return 1.0* np.append(np.arange(1,n+1),np.arange(n+1,0,-1)) / (((n+1)*(n+2))/2 + ((n+1)*(n))/2)
	def expo(n):
		print('not implemented, wma is used')
		return wma(n)
	for i in range(x.size):
		if i<lissage:
			ni = i
		elif x.size-i -1 < lissage:
			ni = x.size-i-1
		else:
			ni = lissage
		if filtre == "lin":
			A = lin(ni)
		elif filtre == "wma":
			A = wma(ni)
		elif filtre == "expo":
			A = expo(ni)
		else:
			A = lin(ni)
		y[i] = np.dot(A,x[i-ni:i+ni+1])

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
	from scipy import fft,ifft
	yfft=fft(x)
	for i in range(x.size):
		if i<x.size/2.0:
			n = 2*i/x.size
		else:
			n = 2*(i-x.size/2.0)/x.size
		if n>lissage:
			yfft[i] = 0.0
	y=ifft(yfft)
	print(y-x)
	return y

def Mean_outliers(x):
	'''from Leys, Christophe, et al. "Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median." Journal of Experimental Social Psychology 49.4 (2013): 764-766.'''
	b = 1.4826
	mx = np.median(x)
	mad = b*np.median(np.abs(x-mx))
	x2 = x[((x<mx+2.0*mad)&(x>mx-2.0*mad))]
	return np.mean(x2)
