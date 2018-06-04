# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 20:06:41 2018

@author: majd
"""

kinds= ['correlation','partial correlation','tangent','covariance']
def freq_dwell(labels,n_states=5,kinds=kinds):
    """
    Computes Group-level frequency of occurence and average dwell time per state. 
    Parameters
    ----------
    labels: array of labels computed from DynamicConnectivityMeasure
    n_states: number of dynamic connectivity states
    kinds: kinds of connecitivty measure 
    
    Returns
    -------
    mean_dwell: Group-level dwell time
    
    counter: Group-level average frequency of occurrence
    
    """
    counter,mean_dwell=dict(),dict()
    for kind in kinds:
        n_subj=labels[kind].shape[0]
        counter[kind]=np.zeros((n_subj,n_states))
        dwell[kind]=np.zeros((n_subj,n_states))
        mean_dwell[kind]=np.empty((n_subj,n_states))
        for i in range(n_subj):
            val_old=1
            for  val in labels[kind][i]:
                for j in range(n_states):
                    if val == j:
                        counter[kind][i,j] = counter[kind][i,j] +1
                    if val != val_old:
                        dwell[kind][i,j]=dwell[kind][i,j]+1
                    val_old=val
                    mean_dwell[kind][i,j]= (counter[kind][i,j]/dwell[kind][i,j])
                    counter[kind][i,j]=100*counter[kind][i,j]/len(labels)
    return mean_dwell,counter 