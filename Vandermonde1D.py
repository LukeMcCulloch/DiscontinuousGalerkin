#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 23:57:29 2019

@author: lukemcculloch
"""
import numpy as np
def Vandermonde1D(N,r):
    """
    % function [V1D] = Vandermonde1D(N,r)
    % Purpose : Initialize the 1D Vandermonde Matrix, V_{ij} = phi_j(r_i);
    #"""
    V1D = np.zeros((len(r),N+1),float)
    for j in range(N+1):
        V1D[:,j] = JacobiP(r[:], 0, 0, j)
    return V1D


if __name__ == """__main__""":
    v = Vandermonde1D(3,[1,1,1])