#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 16:49:41 2018

Convert .npz files to .mat for Vidaurre's HMM toolbox

@author: kosciessa
"""

import numpy as np
import scipy.io as io
import os  
rootPath='/Users/kosciessa/BrainHack/HMM/B_data/10subjects/'
for fn in os.listdir(rootPath):
    dataName=rootPath+fn
    if os.path.isfile(dataName):
        data = np.load(dataName)
        #list(data.keys())
        io.savemat(rootPath+fn[ :-4]+'.mat', mdict={'ts':data['ts'],'id':data['id']})