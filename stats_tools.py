import numpy as np
import la_tools as la
from scipy import stats
from scipy import special
import copy as c

def anova(*args):
	f,p = stats.f_oneway(*args)
	return f,p	
def COV(samples):
	return Covariance_estimate(samples)
def Covariance_estimate(samples):

	try:
		n=samples.shape[1]
		istranspose = True 
	except:
		n=len(samples)
		istranspose = False 
#	print(istranspose)
	if istranspose is True:
		samples2 = samples.T
	else:
		samples2 = c.copy(samples)

	msamples = samples2[0]/n+0.0
	for i in xrange(1,n):
		msamples += samples2[i]/n
#	msamples = msamples/n
	sigma = (1.0/(n-1)) * np.dot(samples2[0][:,np.newaxis],samples2[0][np.newaxis,:])
	for i in xrange(1,n):
#		if i < 5:
#			print(sigma[0:4,0:4])
		sigma += (1.0/(n-1)) * np.dot(samples2[i][:,np.newaxis] - msamples[:,np.newaxis],samples2[i][np.newaxis,:] - msamples[np.newaxis,:])


	return sigma


def Mi2(samples,S=False):
	return Mahalanobis_square(samples,S)
def Mahalanobis_square(samples,S=False):
	try:
		n=samples.shape[1]
		istranspose = True 
	except:
		n=len(samples)
		istranspose = False

	if S is False:
		S = Covariance_estimate(samples)
	if istranspose is True:
		samples2 = samples.T
	else:
		samples2 = c.copy(samples)
	msamples = samples2[0]+0.0
	for i in xrange(1,n):
		msamples += samples2[i]
	msamples = msamples/n
#	print(n)

	Sm1 = np.linalg.pinv(S,rcond=1.0e-12)
	mi2 = np.zeros(n)
	for i in xrange(n):
		mi2[i] = la.mdot(   (  (samples2[i][np.newaxis,:] -msamples[np.newaxis,:] ) , Sm1 , (samples2[i][:,np.newaxis] - msamples[:,np.newaxis])  )   )
	return mi2

def Prob_Mahalanobis(mi):
	data,r = stats.probplot(mi, dist="norm", fit=True,plot=None)
	return data,r

def TMS(samples,S=False):
	return Test_Mardia_skewness(samples,S)
def Test_Mardia_skewness(samples,S=False):
	try:
		n=samples.shape[1]
		d=samples.shape[0]
	except:
		n=len(samples)
		d = samples[0].size
	s = Mardia_skewness(samples,S)
	ddf = d * (d+1.0) * (d+2.0)/6.0
	statistic = n*s/6
	ps = special.gammainc(ddf/2.0,statistic/2.0)
	st = 'skewness: '+ str(s) + ' ddf: ' +str(ddf) +  ' chances: ' + str(100.0*(1.0-2.0*np.abs(0.5-ps)))
	print(st)
	return ps   

def TMK(samples,S=False):
	return Test_Mardia_kurtosis(samples,S)
def Test_Mardia_kurtosis(samples,S=False):
	try:
		n=samples.shape[1]
		d=samples.shape[0]
	except:
		n=len(samples)
		d = samples[0].size

	k = Mardia_kurtosis(samples,S)
	mu = d*(d+2.0)
	sigma2 = 8.0*d*(d+2.0)/n
	sigma2 = sigma2**2
	pk = 0.5 * (   1.0 + special.erf(  (k-mu) /np.sqrt(2.0*sigma2)  )   )
	#pk =  np.exp(  -(k-mu)*(k-mu) /(2.0*sigma2)  )/(np.sqrt(2*np.pi*sigma2))   
	s = 'kurto: '+ str(k) + ', mean th: ' + str(mu) + ', std th: ' + str(np.sqrt(sigma2)) + ' chances: ' +  str(100.0*2.0*(0.5-np.abs(0.5-pk))) + '%'
	print(s)
	return pk


def Mardia_skewness(samples,S=False):
	try:
		n=samples.shape[1]
		istranspose = True 
	except:
		n=len(samples)
		istranspose = False

	if S is False:
		S = Covariance_estimate(samples)
	if istranspose is True:
		samples2 = samples.T
	else:
		samples2 = c.copy(samples)

	Sm1 = np.linalg.pinv(S,rcond=1.0e-12)
	msamples = samples2[0]+0.0
	for i in xrange(1,n):
		msamples += samples2[i]
	msamples = msamples/n

	skew = 0.0
	for i in xrange(0,n):
		for j in xrange(0,n):
			skew += la.mdot(   (  (samples2[i][np.newaxis,:] -msamples[np.newaxis,:] ) , Sm1 , (samples2[j][:,np.newaxis] - msamples[:,np.newaxis])  )   )**3.0
	skew = skew[0,0]/(n*n)
	return skew
def Mardia_kurtosis(samples,mi2=False,S=False):
	if mi2 is False:
		mi2 = Mahalanobis_square(samples,S)
	kurtosis = np.sum(mi2*mi2)/mi2.size
	return kurtosis

def Remove_outliers(mi2,fact=5.0):
	sig = fact*np.std(mi2)
	xi = np.mean(mi2) 
	mi2ct = mi2[mi2<xi+sig]
	mi2c = mi2ct[mi2ct>xi-sig]

	return mi2c

def Sobol(func = False, method = 'salib', MC = False, n_param = False, range_param = False, n_out = False):
	if method == 'salib':
		from SALib.sample import saltelli
		from SALib.analyze import sobol
		problem = {
			'num_vars': 8,
			'bounds': [[0., 1.],
				[0., 1.],
				[0., 1.],
				[0., 1.],
				[0., 1.],
				[0., 1.],
				[0., 1.],
				[0., 1.]]}


def Sobol(func = False, MC = False, QMC = False, n_out=False):
	if func is False:
		print('function needed')
		return
	if MC is False:
		print('MC runs are needed: it is a matrix of the parameters  chosen randomly')
		return
	if QMC is False:
		print('MC runs are needed for the quasi monte carol: it is a matrix of the parameters chosen randomly of the same size than QMC')
		return
	if QMC.shape != MC.shape:
		print('MC and QMC need to be of the same shape')
		return
	if n_out is False:
		n_out = func(MC[0]).size

	n_p = MC.shape[1]
	N = MC.shape[0]

	output_=np.zeros((N,n_out))
	c_out_1=np.zeros((N,n_out,n_p))
	c_out_t=np.zeros((N,n_out,n_p))
	for i in xrange(N):
		input_=MC[i,:]
		output_[i,:] = func(input_)
		for j in xrange(n_p): #pertubation around parameters
			input_=np.r_[QMC[i,0:j],MC[i,j],QMC[i,j+1:]]
			c_out_1[i,:,j]=func(input_) # only one parameter is modified
			input_ = np.r_[MC[i,0:j],QMC[i,j],MC[i,j+1:]]
			c_out_t[i,:,j]=func(input_) # two parameters are modified at the same time

	# monte carlo integrations to estimate integral functions
	f0 = np.sum([output_[x] for x in xrange(N)], axis=0)/N
	# variance of the output
	D = np.sum([output_[x]**2 for x in xrange(N)], axis=0)/N - f0**2
	#D = np.var(output_,axis=0)

	# partial variances
	Dj = np.zeros((n_p,n_out)) # partial variances associated with parameter j
	Dtotj = np.zeros((n_p,n_out)) # total partial variance associated with parameter j
#	Dj = D - np.sum([output_[i] - ], axis = 0)/(2*N)

#	Dtotj = Dtotj/(2*N)
	for i in xrange(N):
		for j in xrange(n_p):
			Dj[j] = Dj[j] + ((output_[i,:] - c_out_1[i,:,j])**2) # start computation of partial variances
			Dtotj[j,:] = Dtotj[j,:] + (output_[i,:] - c_out_t[i,:,j])**2 # total variance due to pj
	Dj = D-Dj/(2*N)
	Dtotj = Dtotj/(2*N)
	#compute sensitivity indices from variances
	print(D)
	Sob_1 = Dj/D # first order
	Sob_t = Dtotj/D # total effect
	print('1st order global Sobol indexes // first order effect :')
	print(Sob_1.flatten())
	print('Total effect Sobol index:')
	print(Sob_t.flatten())
	
	varx = np.zeros(n_p)
	espx = np.zeros(n_p)
	for i in xrange(n_p):
		for j in xrange(N):
			bj = np.copy(MC[j])
			aj = np.copy(QMC[j])
			aij = np.copy(QMC[j])
			aij[i] = np.copy(MC[j][i])
			varx[i] += func(bj) * (func(aij) - func(aj)) 
			espx[i] += (func(aj) - func(aij))**2

	varx = varx/N
	espx = espx/(2*N)

	print(D)
	Sob_1 = varx/D # first order
	Sob_t = espx/D # total effect
	print('1st order global Sobol indexes // first order effect :')
	print(Sob_1.flatten())
	print('Total effect Sobol index:')
	print(Sob_t.flatten())
	return (Sob_1,Sob_t)



















