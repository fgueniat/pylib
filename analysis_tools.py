import numpy as np

def Evolution_operator(K,rtresh = 1.e-15):
    A = np.dot( K[:,1:], np.linalg.pinv(K[:,0:-1],rcond=rtresh))
    return A

def DMDO(K=None,tresh = 0.1, alpha = 0.9,v=False,A=None):

    if A is None:
        if K is None:
            print 'data or Operator needed'
            return -1
        else: A = Evolution_operator(K)
    if v is False:v = np.ones(A.shape[0])
    obs = np.zeros(A.shape[0])
    for i in range(A.shape[0]):
        row = np.zeros(A.shape[0])
        row[np.abs(A[i,:]/v)>tresh] = 1.
        column = np.zeros(A.shape[0])
        column[np.abs(A[:,i]/v[i])>tresh] = 1.
        obs[i] = (alpha * np.sum(row) + (1.-alpha) * np.sum(column))/A.shape[0]

    return obs
 
