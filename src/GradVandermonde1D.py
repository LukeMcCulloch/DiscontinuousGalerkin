#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 18:51:09 2019

@author: luke
"""
import numpy as np
sqrt = np.sqrt

from GradJacobiP import GradJacobiP

def GradVandermonde1D(N,r):
    """        
     function [DVr] = GradVandermonde1D(N,r)
     Purpose : Initialize the gradient of the modal basis (i) at (r) at order N
    """
    DVr = np.zeros( ( len(r),N+1), float)
    
    # Initialize matrix
    for i in range(N+1):
       DVr[:,i] = GradJacobiP( r, 0, 0, i )
    return DVr



if __name__ == """__main__""":
    N=3
    r = np.array([1., 2., 3., 5.])
    vi = GradVandermonde1D(N=N,r=r)
    