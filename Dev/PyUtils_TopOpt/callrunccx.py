#!/usr/bin/env python

from __future__ import print_function
#import os
import numpy as np
from runCCX import runCCX
import readwriteData as rwd


def callrunccx(x,grad,ccxversion,inp,volfrac,rmin,penalty,compfile, plotIterationWiseFlag, scaling):
    
    #plotIterationWiseFlag = False
    #scaling = 1

    ##################################################################
    #Update iteration counter
    callrunccx.counter+=1
    itertop=callrunccx.counter

    #Update filter order
    if (callrunccx.counter) % (callrunccx.q_iter) ==0 and ((callrunccx.q) < (callrunccx.q_max)):
        callrunccx.q*=callrunccx.q_factor


    #Update dnnz values for i=2 onwards
    if callrunccx.counter ==2:
        temp = rwd.maxnnz()
        callrunccx.nnz = int(temp) 

    #Write design variables to density.dat
    rwd.writeCurrentrho(x)

    #Call ccx
    runCCX(ccxversion,inp,penalty,rmin,volfrac,itertop,callrunccx.nnz, callrunccx.q)


    #read rhoPhys and store iterationwise
    if plotIterationWiseFlag:
    	rhoPhys=rwd.getRhoPhys()
    	iterwiseRhofileName = "rho"+str(itertop)+".dat"
    	np.savetxt( iterwiseRhofileName,rhoPhys,fmt='%1.6f') 


    #Get compliance value
    compliance=rwd.getCompl()
    print('Compliance=')
    print(compliance)


    f = open(compfile, "a")
    f.write("\n")
    f.write(str(compliance)) 
    f.close()

    if grad.size > 0:
        grad[:]=rwd.getComplsens()*scaling   
    return compliance*scaling
