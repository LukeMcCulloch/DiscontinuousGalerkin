#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 00:03:10 2019

@author: lukemcculloch
"""
import numpy as np
from scipy.special import gamma , jacobi
sqrt = np.sqrt

def SQ(x):
    return x*x

def JacobiP_scipy(x,alpha,beta,N):
    n = N
    return jacobi(n, alpha, beta, monic=0)

def JacobiP_matlab(x,alpha,beta,N):
    """
    - function [P] = JacobiP(x,alpha,beta,N)
    - Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
    -          (alpha+beta <> -1) at points x for order N and 
                returns P[1:length(xp))]
    - Note   : They are normalized to be orthonormal.
    
    test:
    #
    
    j=3
    r = np.array([1., 2., 3., 5.])
    x=r[:]
    
    alpha=0.
    beta=0.
    N=0
    
    N=2
        
        #
    #"""
    #safety first, a hem:
    aold=0.0
    anew=0.0
    bnew=0.0
    h1=0.0
    gamma0=0.0
    gamma1=0.0
    #
    ab=alpha+beta
    ab1=alpha+beta+1.0
    a1=alpha+1.0
    b1=beta+1.0
    
    
    
    # Turn points into row if needed.
    xp = x
    dims = np.shape(xp)
    if len(dims) ==1:
        dims = (1,dims[0])
    
    if (dims[1]==1):
        xp = xp.conj().transpose()
    
    
    #PL = np.zeros((N+1,np.size(xp)))
    PL = np.zeros((N+1,len(xp)))
    
    # Initial values P_0(x) and P_1(x)
    gamma0 = 2**(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)* \
                    gamma(beta+1)/gamma(alpha+beta+1)
        
    PL[0,:] = 1.0/sqrt(gamma0)
    if N==0 : 
        P=PL.conj().transpose()
        return P
    gamma1 = (alpha+1.)*(beta+1.)/(alpha+beta+3.)*gamma0
    
    PL[1,:] = ( np.asarray((alpha+beta+2)*xp)/2. + (alpha-beta)/2)/sqrt(gamma1)
    if (N==1): 
        P=PL.conj().transpose()
        return P
    
    #TLM current spot:
    # Repeat value in recurrence.
    aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3))
    
    # Forward recurrence using the symmetry of the recurrence.
    for i in range(1,N):
        h1 = 2.*i+alpha+beta
        anew = 2.0/(h1+2.)*sqrt( (i+1.)*(i+1.+alpha+beta)*(i+1.+alpha) * \
                                                   (i+1.+beta)/(h1+1.)/(h1+3.))
        bnew = - (alpha**2 - beta**2)/ h1/(h1+2.) 
        #PL[i+2,:] = 1./anew*( -aold*PL[i,:] + (xp-bnew)*PL[i+1,:] )
        PL[i+1,:] = 1./anew*( -aold*PL[i-1,:] + (xp-bnew)*PL[i,:] )
        aold =anew
        
        
    
    
    P = PL.conj().transpose()
    return P



def JacobiP_cpp(x,alpha,beta,N):
    """
    - function [P] = JacobiP(x,alpha,beta,N)
    - Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
    -          (alpha+beta <> -1) at points x for order N and 
                returns P[1:length(xp))]
    - Note   : They are normalized to be orthonormal.
    
    test:
        #
        j=1
        r = np.asarray([1, 2, 3, 5])
        x=r[:]
        alpha=0
        beta=0
        N=j
        #
    #"""
    #safety first, a hem:
    
    aold   = anew   = bnew = h1 = 0.0
    gamma0 = gamma1 = 0.0
    #
    ab  = alpha+beta
    ab1 = alpha+beta+1.0
    a1  = alpha+1.0
    b1  = beta+1.0
    
    
    x = np.asarray(x)
    Nc = np.size(x)
    
    #    dims = np.shape(x)
    #    if len(dims) ==1:
    #        Nc = (1,dims[0])
        
        
    P = np.zeros((N+1,Nc),float)
    #PL = np.zeros((N+1,np.size(xp)))
    #PL = np.zeros((N+1,len(xp)))
    PL = np.zeros((N+1, Nc),float)
    
    # Initial values P_0(x) and P_1(x)
    gamma0 = ab1**2./(ab1)*gamma(a1)*gamma(b1)/gamma(ab1)
    
    if N==0:
        return 1./sqrt(gamma0)
    else:
        PL[1,:] = 1./sqrt(gamma0)
        
    gamma1 = (a1)*(b1)/(ab+3.0)*gamma0
    prow = ((ab+2.0)*x/2.0 + (alpha-beta)/2.0) / sqrt(gamma1)
    
    if N==1:
        return prow
    else:
        PL[2] = prow
    
    ## Repeat value in recurrence.
    aold = 2.0/(2.0+ab)*sqrt((a1)*(b1)/(ab+3.0))
    
    ## Forward recurrence using the symmetry of the recurrence.
    for i in range(1,N):
        h1 = 2.0*i+ab
        anew = 2.0/(h1+2.0)*sqrt((i+1)*(i+ab1)*(i+a1)*(i+b1)/(h1+1.0)/(h1+3.0))
        bnew = - (SQ(alpha)-SQ(beta))/h1/(h1+2.0)
        x_bnew = x-bnew
        PL[i+2,:] = 1.0/anew*( -aold*PL[i] + x_bnew * PL[i+1] )
        aold =anew 
    
    return P