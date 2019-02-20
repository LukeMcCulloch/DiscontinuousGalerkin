#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 00:03:10 2019

@author: lukemcculloch
"""
import numpy as np
from scipy.special import gamma #, factorial
sqrt = np.sqrt

def JacobiP(x,alpha,beta,N):
    """
    - function [P] = JacobiP(x,alpha,beta,N)
    - Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
    -          (alpha+beta <> -1) at points x for order N and 
                returns P[1:length(xp))]
    - Note   : They are normalized to be orthonormal.
    #"""
    #safety first, a hem:
    aold=0.0
    anew=0.0
    bnew=0.0
    h1=0.0
    gamma0=0.0
    gamma1=0.0
    # Turn points into row if needed.
    xp = x
    dims = np.shape(xp)
    if (dims(1)==1):
        xp = xp.conj().transpose()
    
    
    PL = np.zeros(N+1,np.size(xp))
    
    # Initial values P_0(x) and P_1(x)
    gamma0 = 2**(alpha+beta+1)/(alpha+beta+1)*gamma[alpha+1]* \
        gamma[beta+1]/gamma[alpha+beta+1]
        
    PL[1,:] = 1.0/sqrt(gamma0)
    if (N==0): 
        P=PL.conj().transpose()
        return
    gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0
    PL[2,:] = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/sqrt(gamma1)
    if (N==1): 
        P=PL[N+1,:].conj().transpose()
        return
    
    #TLM current spot:
    # Repeat value in recurrence.
    aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3))
    
    # Forward recurrence using the symmetry of the recurrence.
    for i in range(N-1):
        h1 = 2*i+alpha+beta
        anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)* \
          (i+1+beta)/(h1+1)/(h1+3))
        bnew = - (alpha**2-beta**2)/h1/(h1+2)
        aold =anew
    
    P = PL[N+1,:].conj().transpose()
    return P