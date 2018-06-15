#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 20:41:15 2018

Load HMM output that has been exported from MATLAB

The following parameters have been used:
    12 states estimated across subjects
    individual state occupancies etc.
    unique full covariance matrix to quench interindividual differences of mean connectivity matrices
    For more details see https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#

@author: kosciessa
"""

import scipy.io
file='/Users/kosciessa/BrainHack/HMM/B_data/A_HMMoutput.mat'
mat = scipy.io.loadmat(file)