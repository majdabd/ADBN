# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 19:12:13 2018

@author: majd
"""
import nilearn.datasets as dataset
from nilearn.input_data import NiftiLabelsMasker , NiftiMapsMasker , NiftiSpheresMasker
from nilearn.datasets import load_mni152_template,load_mni152_brain_mask
brainmask = load_mni152_brain_mask()
mem = Memory('nilearn_cache')

def signal_extract(data,atlas,t_r=2.2,masker_type='Spheres',saveas='file'):
    
    """
    Extracts BOLD time-series from regions of interest
    
    Parameters
    ----------
    data: Filenames of subjects. 
    
    atlas: regions or coordinates to extract signals from.
    
    masker_type : Type of masker used to extract BOLD signals . types are : 'Spheres','Maps','Labels'
    
    saveas : Destination to save and load output (.npz)
    
    Returns
    ---------
    subject_ts : array-like , 2-D (n_subjects,n_regions)
                 Array of BOLD time-series 
    """
    subjects_ts=[]
    
    if os.path.exists(saveas):
        
        subjects_ts=np.load(saveas)['arr_0']
        
    else:
        
        if masker_type== 'Spheres':
            masker = NiftiSpheresMasker(
                            seeds=atlas, smoothing_fwhm=6, radius=4 ,mask_img=brainmask,
                            detrend=False, standardize=True, low_pass=0.1, high_pass=0.01, t_r=t_r)
        elif masker_type == 'Maps':
            masker = NiftiMapsMasker(maps_img=atlas,mask_img=brainmask,standardize=True,
                                 high_pass=0.01,low_pass=0.1,detrend=False,t_r=t_r,
                                 memory_level=2,smoothing_fwhm=5,resampling_target='data',
                                 memory=mem,verbose=5)
        elif masker_type == 'Labels':
            masker = NiftiLabelsMasker(labels_img=atlas,mask_img=brainmask,standardize=True,
                                 high_pass=0.01,low_pass=0.1,detrend=False,t_r=t_r,
                                 memory_level=2,smoothing_fwhm=5,resampling_target='data',
                                 memory=mem,verbose=5)
        else:
            raise ValueError("Please provide masker type")
            
        for func_file in data:
            time_series = masker.fit_transform(func_file)
            subjects_ts.append(time_series)
            np.savez(saveas,subjects_ts)
            
    return subjects_ts
