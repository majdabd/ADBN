# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 20:06:41 2018

@author: majd
"""

kinds= ['correlation','partial correlation','tangent','covariance']
mods=['correlation']
def freq_dwell(labels,n_states=5,mods=mods):
    counter,trans,dwell,freq,grp_freq,grp_dwell=={},{},{},{},{},{}
    for mod in mods:
        n_subj=labels[mod].shape[0]
        counter[mod]=np.zeros((n_subj,n_states))
        trans[mod]=np.zeros((n_subj,n_states))
        dwell[mod]=np.zeros((n_subj,n_states))
        freq[mod]=np.zeros((n_subj,n_states))
        for i in range(n_subj):
            val_old=0
            for  val in labels[mod][i]:
                for j in np.arange(0,n_states):
                    if val == j:
                        counter[mod][i,j] = counter[mod][i,j] +1          

                    if val != val_old and val==j:
                        trans[mod][i,j] += 1
                val_old=val          
        dwell[mod]= counter[mod]/trans[mod]
        freq[mod]=counter[mod]*100/labels[mod][i].shape[0]
        grp_freq[mod]=freq[mod].mean(axis=0)
        grp_dwell[mod]=dwell[mod].mean(axis=0)
    return dwell,freq,grp_dwell,grp_freq
