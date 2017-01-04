from sklearn.neighbors import NearestNeighbors
import numpy as np

def Find_neighbors(x,X,k=20,nbrs = False):
	if X.shape[1] ==0:
		NN = np.array([])
	else:
		if X.shape[1]<k:
			NN = np.arange(X.shape[1])
		else:	
			if nbrs is False:
				nbrs = index(X,k)
			dist,NN = nbrs.kneighbors(x.reshape(1, -1))
	return NN.flatten()

def index(X,k=20):
	if X.shape[0]>1000:
		nbrs = NearestNeighbors(n_neighbors=k, algorithm='auto',n_jobs=-1).fit(X.T)
	else:
		nbrs = NearestNeighbors(n_neighbors=k, algorithm='auto').fit(X.T)
	return nbrs

def KNN(x,X,k=20,nbrs=False):
	return Find_neighbors(x,X,k,nbrs)


