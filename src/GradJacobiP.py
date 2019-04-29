#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 18:51:38 2019

@author: luke
"""
import numpy as np
sqrt = np.sqrt

from JacobiP import JacobiP_scipy, JacobiP_matlab, JacobiP_cpp
JacobiP = JacobiP_matlab

def GradJacobiP(r, alpha, beta, N):
    """function [dP] = GradJacobiP(z, alpha, beta, N);
        Purpose: Evaluate the derivative of the orthonormal Jacobi
  	   polynomial of type (alpha,beta)>-1, at points x
            for order N and returns dP[1:length(xp))]
    """
    
    dP = np.zeros( len(r), float)
    #dP = np.zeros( (np.size(r), 1), float)
    #dP = np.zeros( (np.shape(r)[1], 1), float)
    if(N == 0):
        dP[:] = 0.0
    else:
        dP[:] =  sqrt(N*(N+alpha+beta+1)) * JacobiP(r,alpha+1,beta+1, N-1)[:,N-1] 
    return dP