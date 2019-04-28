#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 11:42:40 2019

@author: lukemcculloch
"""
class SolveProp(object):
    #Parameters
    gamma = 1.4; CFL = 1.0; time = 0;
    
    # Prepare for adaptive time stepping
    mindx = min(x(2,:)-x(1,:))

    # Limit initial solution
    rho =SlopeLimitN(rho)
    rhou=SlopeLimitN(rhou) 
    Ener=SlopeLimitN(Ener)


def Euler1D(rho,rhou,Ener,FinalTime):
    
    return