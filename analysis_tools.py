import numpy as np
#from scipy.cluster.vq import kmeans2 as kmeans
from scipy.cluster.vq import whiten as precond
from scipy.cluster.vq import kmeans2
from sklearn.cluster import KMeans as kmeans
from sklearn.metrics import pairwise_distances_argmin_min as pairwise_dist
import la_tools as lat
from scipy.signal import detrend

def evolution_operator(K,rtresh = 1.e-15,ismean=True):
    if ismean is False:
        K2 = detrend(K,axis=1)
        A = np.dot( K2[:,1:], np.linalg.pinv(K2[:,0:-1],rcond=rtresh))
    if ismean is True:
        A = np.dot( K[:,1:], np.linalg.pinv(K[:,0:-1],rcond=rtresh))
    return A

def compute_doc(K=None,tresh = 0.1, alpha = 0.9,v=False,A=None,normalize = False):
    ''' 
    Compute the observability criterion out of data or of the flow map operator
    '''
    if A is None:
        if K is None:
            print('data or flow map Operator needed')
            return -1
        else: A = evolution_operator(K)
    if v is True:
        if K is not None:
            v = np.mean(K,axis=1)
        else:
            print('provide K for average renormalization')
            v = False

    if v is False:v = np.ones(A.shape[0])
    v[np.abs(v)<1.e-8] = 1.    
    obs = np.zeros(A.shape[0])
    for i in range(A.shape[0]):
        row = np.zeros(A.shape[0])
        row[np.abs(A[i,:]/v)>tresh] = 1.
        column = np.zeros(A.shape[0])
        column[np.abs(A[:,i]/v[i])>tresh] = 1.
        obs[i] = (alpha * np.sum(row) + (1.-alpha) * np.sum(column))/A.shape[0]
    if normalize is True:
        obs = obs - np.min(obs)        
        obs = obs/np.max(obs)
        indzero = np.where(obs==0)[0]
        do =  np.min(obs[np.nonzero(obs)])
        fact = indzero.size/obs.size
        obs = obs - do * (1. - fact)
        obs[obs<0.] = 0.
        obs = obs/np.max(obs)
    return obs
 



def cluster_doc(doc=None,grid_x=None,grid_y=None,nx=None,ny=None,nclust=10,renorm=True,recluster=False,ncluster_recluster = False):
    '''
    Function for identifying main observervable areas.
    doc: observability estimate. Call DMDO.
    '''
    if doc is None:
        print('please provide the criterion by calling function compute_doc')
        return -1
    N=doc.size
    n=1
    if grid_x is not None:n=n+1
    if grid_y is not None:n=n+1
    if n ==1: 
        print('please provide grid')
        return -1
    if grid_x is None and n==2:
        return cluster_obs(doc=doc,grid_x=grid_y,grid_y=None,nx=nx,ny=ny,nclust=nclust,renorm=renorm,recluster = recluster,ncluster_recluster=ncluster_recluster)
    if n==3:
        if ny is None:
            try:
                ny = grid_y.shape(1)
            except:
                print('please provide either 2D grids or shapes of the grid')
                return -1
        if nx is None:
            try:
                nx = grid_x.shape(0)
            except:
                print('please provide either 2D grids or shapes of the grid')
                return -1
        if nx * ny != grid_x.size:
            print('errors in shape')
            return -1
    if n == 2:
        if nx is None:
            nx = grid_x.shape(0)
        if nx != grid_x.size:
            print('errors in shape')
            return -1
        ny=1

    if recluster is True:
        if ncluster_recluster is False:
            ncluster_doc= 1+int(nclust/2.)
        else:
            ncluster_doc = ncluster_recluster
        cent,labels = kmeans2(doc,ncluster_doc)
        perm = np.argsort(cent)
        doc_t = np.array([ perm[l] for l in labels ])/(ncluster_doc-1.)
    else:
        doc_t = np.copy(doc)    
    
    data = np.zeros((N,n))
    if np.linalg.norm(doc)>1.e-8:
        data[:,0] = doc.flatten()
        #1.*(doc-np.min(doc))/(np.max(doc)-np.min(doc))
    else:
        print('doc is null')
        return -1
    for k in range(nx*ny):
        ix,iy = np.unravel_index(k,(nx,ny))
        if n==3: x = grid_x[ix,iy]
        else: x = grid[ix]
        if n==3: y = grid_y[ix,iy]
        data[k,1] = x
        if n==3: 
            data[k+nx*ny,1] = x
            data[k,2] = y
            data[k+nx*ny,2] = y
    if renorm is True:
        #clusters,keys = kmeans(precond(data),nclust)
        computed_clusters = kmeans(n_clusters=nclust).fit(precond(data))
    else:
        #clusters,keys = kmeans(data,nclust)
        computed_clusters = kmeans(n_clusters=nclust).fit(data)
    
    isperm=False
    doc_c = np.zeros(nclust)
    if isperm is True:
        for iclust in range(nclust):
            ind_clust = np.where(computed_clusters.labels_ == iclust)
            doc_c[iclust] = np.mean(doc_t[ind_clust[0]])  
        perm = np.argsort(doc_c)
    else:
        perm = np.arange(nclust)
    keys = np.array(
            [perm[key] for key in computed_clusters.labels_]
            )
    
    clusters = computed_clusters.cluster_centers_
    if renorm is False:
        clusters_closest_t, _ = pairwise_dist(clusters, data)
    else:
        clusters_closest_t, _ = pairwise_dist(clusters, precond(data))
    
    clusters_closest = np.array(
                [clusters_closest_t[perm[p]] for p in perm]
                ).astype(int)
    clusters_max_t = np.zeros(nclust)
    docs_t = np.zeros(nclust)
    for key in range(nclust):
            ind_key, = np.where(computed_clusters.labels_ == key)
            doc_key = doc[ind_key]
            ind_max = np.argmax(doc_key)
            docs_t[key] = doc_key[ind_max]
            clusters_max_t[key] = ind_key[ind_max]
    print(docs_t)
    docs = np.array([ docs_t[perm[key]] for key in range(nclust)])
    clusters_max = np.array([ clusters_max_t[perm[key]] for key in range(nclust)]).astype(int)
    
    return clusters_max,docs,keys,clusters_closest


