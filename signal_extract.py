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

def signal_extract(func_data=None,confounds=None,atlas_img=None,masker_type='Spheres',smoothing_fwhm=6,high_pass=0.01,low_pass=0.1,t_r=2.2,detrend=False,saveas='file'):
    
    """
    Extracts BOLD time-series from regions of interest
    
    Parameters
    ----------
    func_data: functional images ( Default= None ) 
    
    confounds: Confounds file used to clean signals ( Default= None )
    
    atlas_img: regions or coordinates to extract signals from ( Default= None )
    
    masker_type : Type of masker used to extract BOLD signals . types are : 'Spheres','Maps','Labels'
    
    smoothing_fwhm : Smoothing width applied to signals in mm ( Default= 6 mm )
    
    high_pass, low_pass: Bandpass-Filtering ( Default= 0.01-0.1 Hz )
    
    detrend: Detrending signals ( Default= False )
    
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
        
        if 
        
        if masker_type== 'Spheres':
            masker = NiftiSpheresMasker(
                            seeds=atlas_img, smoothing_fwhm=smoothing_fwhm, radius=4 ,mask_img=brainmask,
                            detrend=False, standardize=True, low_pass=low_pass, high_pass=high_pass, t_r=t_r
            )
        elif masker_type == 'Maps':
            masker = NiftiMapsMasker(
                                    maps_img=atlas_img,mask_img=brainmask,standardize=True,
                                    low_pass=low_pass, high_pass=high_pass, t_r=t_r,
                                    memory_level=2,smoothing_fwhm=smoothing_fwhm,resampling_target='data',
                                    memory=mem,verbose=5
            )
        elif masker_type == 'Labels':
            masker = NiftiLabelsMasker(
                                 labels_img=atlas_img,mask_img=brainmask,standardize=True,
                                 high_pass=high_pass,low_pass=low_pass,detrend=False,t_r=t_r,
                                 memory_level=2,smoothing_fwhm=smoothing_fwhm,resampling_target='data',
                                 memory=mem,verbose=5
            )
            
        else:
            raise ValueError("Please provide masker type")
        
        if confounds is not None:    
            for func_file, confound_file in zip(func_data,confounds):
                time_series = masker.fit_transform(func_file,confounds=confound_file)
                subjects_ts.append(time_series)
                np.savez(saveas,subjects_ts)
        else:
          for func_file in data:
            time_series = masker.fit_transform(func_file)
            subjects_ts.append(time_series)
            np.savez(saveas,subjects_ts)   
            
    return subjects_ts
