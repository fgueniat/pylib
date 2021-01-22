from sklearn.neighbors import NearestNeighbors,KNeighborsClassifier
import stats_tools as st
import peaks as pks
import numpy as np
from scipy.cluster.vq import kmeans2 as kmeans

####
# kmeans !
####

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

def find_peaks(y_axis, x_axis = None, lookahead = 200, delta=0):
    r1,r2 = pks.peakdetect(y_axis, x_axis = None, lookahead = 200, delta=0)
    return np.array(r1),np.array(r2)


class result():
    def __init__(self,centroids=None,classifier = None,nc = None,training_labels = None):
        self.centroids = centroids
        self.classifier = classifier
        self.nc = nc
        self.training_labels = training_labels
        self.entries_labels = None
    
    def reinit(self,res2):
        self.centroids = res2.centroids
        self.classifier = res2.classifier
        self.nc = res2.nc
        self.training_labels = res2.training_labels

    def cluster(self,input_):
        try:
            if len(input_.shape) == 1:
                input_ = input_.reshape(-1,1)
            return self.classifier.predict(input_)[0]
        except:
            print('classifier not set or try input.reshape(-1,1)')
            return -1
    def id_newlabels(self,entries):
        self.entries_label = [self.cluster(input_) for input_ in entries]

def kmeans_classifier(training = None, entries = None, centroids = None, classifier=None, nc=10, res=None):
    if res is None:
        res = result(centroids=centroids,classifier=classifier,nc=nc)
    
    if res.classifier is None:
        if training is None:
            print('data needed to build the classifier')
            return -1
        else:
            try:
                res.centroids,res.training_labels = kmeans(training,k=nc,missing='raise')
            except:
                res.centroids,res.training_labels = kmeans(training,k=nc)
            res.classifier = KNeighborsClassifier()
            res.classifier.fit(training,res.training_labels)
        #
    #
    if entries is not None:
        res.id_newlabels(enties)
    return res
    





