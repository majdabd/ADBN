# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 19:14:50 2018

@author: majd
"""
import os
import nilearn.connectome
import numpy as np
from scipy import signal
import scipy.linalg 
from scipy.linalg import eigh
from sklearn.cluster import KMeans,SpectralClustering,DBSCAN
from sklearn.covariance import EmpiricalCovariance,LedoitWolf,ShrunkCovariance,OAS
from pyriemann.clustering import Kmeans
from pyriemann.utils.base import sqrtm, invsqrtm, logm, expm
import statsmodels.stats.api as sms
from bct.algorithms.modularity import modularity_und
from bct.algorithms.modularity import community_louvain
from bct import grid_communities
def sliding_window_1d(a,ws=60,ss=1):
    """
    Parameters
    ----------
    a:  array , 1-D (n_volumes)
    ws: window width
    ss: window step size, in samples. If not provided, window and step size are equal.
    
    Returns
    ----------
    out : array-like , 2-D (n_windows,ws)
          Vector with all sliding windows
    """ 
    if None is ss:
    # no step size was provided. Return non-overlapping windows
        ss = ws   
    # calculate the number of windows to return, ignoring leftover samples, and
    # allocate memory to contain the samples
    valid = len(a) - ws
    
    nw = abs((valid) // ss)
    
    out = np.ndarray((nw,ws),dtype = a.dtype)
    
    for i in range(nw):
        
        # "slide" the window along the samples
        start = i * ss
        stop = start + ws
        out[i] = a[start : stop]
    
    return out

def taper(window,ws=60,wtype='tukey'):
    """
    Parameters
    ----------
    window : array-like , array of sliding window indices 
    ws: window width
    wtype: type of tapering function ( default : Tukey )
    
    Returns
    -------
    out :
    
    """
    if wtype == 'tukey':
        taper= scipy.signal.tukey(ws,alpha=0.5,sym=True)
      
    out=window*taper
    
    return out

def upper_tri(matrix):
        mat_triu=[]
        for i,x in enumerate(matrix):
                cur_triu = (np.triu(x,1))
                mat_triu.append(cur_triu[np.triu_indices_from(cur_triu,1)])
        return mat_triu
    
def unupper_tri(X):
    n_features=int((1+np.sqrt(1+8*X.shape[1]))//2)
    centroids = np.zeros((X.shape[0],n_features,n_features))
    for i,center in enumerate(X):
            curmat=np.zeros((n_features,n_features))
            curmat[np.triu_indices_from(curmat,1)]=center
            centroids[i] = curmat+curmat.T
    return centroids
    
        
def is_spd(A, tol=1e-8):
    """
    Parameters
    ----------
    A : Symmetrtic Square Matrix 
    
    Returns
    ---------
    Boolean True or False
    """
    E,V = scipy.linalg.eigh(A)
    return np.all(E > -tol),E

def untangent_space(T, Cref):
    """Project a set of Tangent space vectors in the manifold according to the given reference point Cref
    
    Parameters
    ----------
    T:    {array-like} ,The Tangent space , shape= ( NWindows X Mfeatures(*(Mfeatures-1)/2))
    
    Cref: {array-like} ,The reference covariance matrix (Mfeatures X Mfeatures )
    
    Returns
    ----------
    
    covmats: {array-like} ,(Mfeatures X Mfeatures)  SPD Matrice

    """
    Nt, Nd = T.shape
    Ne = (1+np.sqrt(1+8*Nd))/2
    C12 = sqrtm(Cref)
    idx = np.triu_indices_from(Cref,0)
    covmats = np.zeros((Nt, Ne, Ne))
    covmats[:, idx[0], idx[1]] = T
    for i in range(Nt):
        covmats[i] = np.diag(np.diag(covmats[i])) + np.triu(
            covmats[i], 0) + np.triu(covmats[i], 0).T 
        covmats[i] = expm(covmats[i])
        covmats[i] = np.dot(np.dot(C12, covmats[i]), C12)
    return covmats

def conn_highvariance(allcovdata):

    """ 
    -Identify windows with high variance in connectivity for each subject, 
    -calculate the average connectitivy (average of all edges)
    -Define a 95% confidence interval on this average 
    -Select data points outside (higher values)
    
    Parameters
    ----------
    allcovdata:{array-like} , Connectivity matrices of all subjects , 
                shape =(Subjects X NWindows X Mfeatures X Mfeatures )
    
    Returns
    ----------
    
    mtd_allsubj_highvar:{array-like} , All windows of High Variance , shape= (Windows X Mfeatures X Mfeatures)

    """
    mtd_allsubj_highvar = []
    
    # High variance windows for each subject 
    var_mtd_allsubj =[]
    
    for curcov in allcovdata:

        # calculate variance of MTD intra subject
        var_mtd_allsubj.append([a.mean() for a in curcov])

    # Extract points with high variance ( > 95 % confidence interval )
    for cur_i,curvarmtd in enumerate(var_mtd_allsubj):

        a = sms.DescrStatsW(curvarmtd)
        _,high= a.tconfint_mean()
        ind_highvar = np.argwhere(curvarmtd>high)

        # select the covdata for these points only
        curcov = allcovdata[cur_i]
        mtd_allsubj_highvar.append(curcov[ind_highvar])
    return np.vstack(mtd_allsubj_highvar)

def dfc_slid_window(X,ws,ss=1): 
        """
        Computes Sliding-window time-series per subject per region.
        Applies a Tukey window for tapering 
        
        Parameters
        ----------
        X: {array-like}, shape = (n_subjects , n_volumes , n_regions)
           resting state BOLD time-series for subjects
        ws : Window size 
        ss: Sliding step size (default=1)
        
        Returns
        -------
        slwin_ts : Array-like (n_subjects,n_windows,ws,n_regions)
        """
        nsubj=X.shape[0] # number of subjects
        nvolm=X.shape[1] # number of volumes
        nfeat=X.shape[2] # number of brain regions   
        slwin_ts=np.ndarray((nsubj,np.int16(np.ceil((nvolm-ws)//ss)),ws,nfeat))
        for idx,s in enumerate(X):
            fulltimewin = np.arange(nvolm,dtype='int32')
            swins= sliding_window_1d(a=fulltimewin,ws=ws,ss=ss)
            n_slwin = swins.shape[0] #number of sliding windows
            slwin_ts[idx]=np.empty((n_slwin,ws,nfeat))
            for n, curwin in enumerate(swins):
                cur_ts = s[curwin,:]
                slwin_ts[idx][n]=np.ndarray((ws,nfeat))
                for i in range(nfeat):
                    slwin_ts[idx][n][:,i]= taper(cur_ts[:,i],ws)
        return slwin_ts,slwin_ts.shape[1]


class DynamicConnectivityKmeans():
    
    """ Dynamic Connectivity Estimator using two passes of Kmeans on FC matrices
    Parameters
    -----------
    cov_estimator : Method of Estimating covariance/connectivity ( default= EmpiricalCovariance)
    n_states: Number of Dynamic Connectivity States
    ws: sliding window width
    ss: sliding window step size 
    kinds : list of connectivity measures
    
    Attributes
    -----------
    connectivity: {array-like}, Windowed Connectivity Measures for Subjects , 
                  shape=(n_subjects X n_windows X n_features X n_features)
    states : {array-like} ,Dynamic States , shape=(n_states , n_features , n_features)
    labels : {array-like} , Kmeans Labels, shape=(1 , n_windows)
    inertia : Kmeans inertia value
    mean_: {array_like}, Geometric mean for the tangent kind , shape=(n_features , n_features)
    
    """
    def __init__(self,cov_estimator=EmpiricalCovariance(store_precision=True, 
                 assume_centered=False),n_states=2,ws=60,ss=1,kinds=['correlation'],
                 saveas='file'):
        
        self.cov_estimator=cov_estimator
        self.n_states=n_states
        self.ws=ws
        self.ss=ss
        self.kinds=kinds        
        self.saveas=saveas
        
    
        
    def fit(self,X):
        
        """
        Parameters
        -----------
        X: {array-like}, shape = (n_subjects , n_volumes , n_features)
           resting state BOLD time-series for subjects
        
        Returns
        ----------
        self : object 
        
        """    
        if not os.path.exists(self.saveas):
            
            #Checking if input Time series Matrix is 2D for each subject
            subjects_types = [type(s) for s in X]
            if set(subjects_types) != set([np.ndarray]):
                raise ValueError("Each subject must be 2D numpy.ndarray.\n You "
                                 "provided {0}".format(str(subjects_types)))
                
            #Initializing   
            self.connectivity,self.states,self.labels,self.inertia=dict(),dict(),dict(),dict()
            connectivity=dict()
            nsubj=X.shape[0] # number of subjects
            nvolm=X.shape[1] # number of volumes
            nfeat=X.shape[2] # number of brain regions  
            ws=self.ws
            ss=self.ss
            
            #Estimating Dynamic Connectivity
            for kind in self.kinds:
                estimator = nilearn.connectome.ConnectivityMeasure(kind=kind,cov_estimator=self.cov_estimator,
                                                                   vectorize=True,discard_diagonal=True)
                
                slwin_ts,n_slwin=dfc_slid_window(X,ws=ws,ss=ss) # Windowed BOLD time-series
                connectivity[kind] = np.ndarray((nsubj*n_slwin,(nfeat*(nfeat-1))//2))     
                connectivity[kind] = estimator.fit_transform(slwin_ts.reshape(nsubj*n_slwin,ws,nfeat))
                
                #Extracting matrices with high-variance
                highvar_mats=connectivity[kind].reshape(nsubj,n_slwin,(nfeat*(nfeat-1))//2)
                highvar_mats=np.vstack(conn_highvariance(highvar_mats))
                
                
                # 2-level KMeans clustering
                clustering_lvl_1=KMeans(n_clusters=self.n_states, init='k-means++', n_init=150, max_iter=300, 
                                       precompute_distances='auto', verbose=0, random_state=None, 
                                       copy_x=True, algorithm='auto').fit(highvar_mats) # 1st level
                
                clustcenters=clustering_lvl_1.cluster_centers_ # 1st level cluster centers
                
                clustering_lvl_2=KMeans(n_clusters=self.n_states, init=clustcenters, n_init=150, max_iter=300, 
                                       precompute_distances='auto', verbose=0, random_state=None, 
                                       copy_x=True, algorithm='auto').fit(connectivity[kind]) # 2nd level
                
                
                #Finalizing
                self.connectivity[kind]= unupper_tri(connectivity[kind])
                self.states[kind] =unupper_tri(clustering_lvl_2.cluster_centers_)
                self.labels[kind]= clustering_lvl_2.labels_.reshape(nsubj,n_slwin)
                self.inertia[kind]=clustering_lvl_2.inertia_
                
                #Projecting from Tangent space to SPD for the tangent kind 
                if kind=='tangent':
                    self.mean_=unupper_tri(estimator.mean_)
                    self.states[kind]=untangent_space(self.states[kind],Cref=self.mean_) 
                    
            np.savez(self.saveas,self)    
        
        else:
            data=np.load(self.saveas)
            self=data['arr_0'].flatten()[0] 
                        
        return self
