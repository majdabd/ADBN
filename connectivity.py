# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 19:14:20 2018

@author: majd
"""

from sklearn.covariance import LedoitWolf , OAS , GraphLassoCV , EmpiricalCovariance
kinds=['tangent','correlation','partial correlation','covariance']

def connectivity(subjects_ts,kinds=kinds,saveas='file'):
    
    """
    Estimates Functional Connectivity using several estimation models 
    Parameters
    ----------
    subjects_ts: array-like , 2-D (n_subjects,n_regions)
                 Array of BOLD time-series  
    
    kinds: list of kinds of connectivity measure to be computed . kinds include : 
        ' correlation ' , ' partial correlation', ' tangent' , 'covariance' . 
                                                
    
    saveas : Destination to save and load output (.npz)
    
    Returns
    ---------
    mean_connectivity_matrix: dictionary ,  {'kind' : (n_regions,n_regions)} 
                              Group-level functional connectivity matrix
    individual_connectivity_matrix: dictionary , {'kind' : (n_subjects,n_regions,n_regions)}
                              Subject-level functional connectivity matrices
                 
    """
    
    individual_connectivity_matrices = dict()
    
    mean_connectivity_matrix = dict()
    
    if os.path.exists(saveas):
            data=np.load(saveas)
            individual_connectivity_matrices=data['arr_0'].flatten()[0] 
            mean_connectivity_matrix=data['arr_1'].flatten()[0]
    else:

        for kind in kinds:
            
            # Computing individual functional connectivity
            
            conn_measure = ConnectivityMeasure(cov_estimator=LedoitWolf(
                                               assume_centered=True, 
                                               store_precision=True), 
                                               kind=kind, 
                                               vectorize=False, discard_diagonal=False)
            
            individual_connectivity_matrices[kind]= conn_measure.fit_transform(subjects_ts)
            
            # Computing group functional connectivity
            
            if kind == 'tangent':
                mean_connectivity_matrix[kind] =  conn_measure.mean_
            else:
                mean_connectivity_matrix[kind] = \
                individual_connectivity_matrices[kind].mean(axis=0)
            np.savez(saveas,individual_connectivity_matrices,mean_connectivity_matrix)
            
    return mean_connectivity_matrix,individual_connectivity_matrices