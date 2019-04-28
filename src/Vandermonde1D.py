#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 23:57:29 2019

@author: lukemcculloch
"""
import numpy as np

from JacobiP import JacobiP_scipy,JacobiP_matlab, JacobiP_cpp

JacobiP_ = JacobiP_matlab

#JacobiP_ = JacobiP_cpp

def Vandermonde1D(N,r):
    """
    - function [v1d] = vandermonde1d(n,r)
    - purpose : initialize the 
        1d vandermonde matrix, v_{ij} = phi_j(r_i);
    #"""
    v1d = np.zeros((len(r),N+1),float)
    for j in range(N+1):
        v1d[:,j] = JacobiP_(x=r[:], 
                               alpha=0, 
                               beta=0, 
                               N=j)[:,j]
    return v1d


def Vandermonde1D_scipy(N,r):
    """
    - function [V1D] = Vandermonde1D(N,r)
    - Purpose : Initialize the 
        1D Vandermonde Matrix, V_{ij} = phi_j(r_i);
    #"""
    return np.vander(r, N)


if __name__ == """__main__""":
    N=3
    r = np.array([1., 2., 3., 5.])
    vi = Vandermonde1D(N=N,r=r)
    
    v = Vandermonde1D_scipy(N, r)