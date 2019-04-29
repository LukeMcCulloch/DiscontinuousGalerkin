#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 21:30:37 2019

@author: luke
"""
import numpy as np
from GradVandermonde1D import GradVandermonde1D


def Dmatrix1D(N,r,V):
    """
        Dr = Dmatrix1D(N,r,V)
        Purpose : Initialize the (r) 
                    differentiation matrices on the interval,
                    evaluated at (r) at order N
    """
    
    Vr = GradVandermonde1D(N, r);
    Dr = np.linalg.inv(V)*Vr
    return Dr



if __name__ == """__main__""":
    N=3
    r = np.array([1., 2., 3., 5.])
    vi = GradVandermonde1D(N=N,r=r)
    
    