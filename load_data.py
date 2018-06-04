# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 19:08:41 2018

@author: majd
"""

def load_data(subj_list,location='/Users/'):
    """
    Loads resting state functional Images and corresponding confounds
    
    Parameters
    ----------
    subj_list: list of subjects' names
    
    Returns
    -------
    filename: list of resting-state data
    confounds: list of confounds 
    """
    filename=[nib.load('%s/%s_task-rest_preproc.nii.gz'%(location,name)) 
                for name in subj_list]
    confounds=['%s/%s_confounds.csv'%(location,name)]
    return filename,confounds