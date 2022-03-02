#!/usr/bin/env python

from __future__ import print_function
#import os
import numpy as np
import subprocess
import readwriteData as rwd


def callrunccx(x,grad,ccxversion,inp,volfrac,rmin,penalty,compfile,  plotIterationWiseFlag, scaling, passivefilename):
    """
    Call ccx and calculate compliance
    """   

    #Update iteration counter
    callrunccx.counter+=1
    itertop=callrunccx.counter

    #Update filter order
    if (callrunccx.counter) % (callrunccx.q_iter) ==0 and ((callrunccx.q) < (callrunccx.q_max)):
        callrunccx.q*=callrunccx.q_factor

    #Update beta for Heaviside
    callrunccx.beta = rwd.update_beta(callrunccx.counter, callrunccx.beta_iter_init, callrunccx.beta_iter, callrunccx.beta, callrunccx.beta_factor, callrunccx.betamax)


    #Update dnnz values for i=2 onwards
    if callrunccx.counter ==2:
        temp = rwd.maxnnz()
        callrunccx.nnz = int(temp) 

    #Write design variables to density.dat
    rwd.writeCurrentrho(x)

    #Call ccx
    runCCX(ccxversion,inp,penalty,rmin,volfrac,itertop,callrunccx.nnz, callrunccx.q, callrunccx.beta, callrunccx.eta)


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
        g=rwd.getComplsens()

        if isinstance(passivefilename, str):
            passiveEl = rwd.readElemList(passivefilename)
            g = rwd.setElemVal(g, passiveEl,0* np.ones(np.size(passiveEl))) 

        grad[:]=g*scaling   
    return compliance*scaling


def getVolconstraint(x, grad,ccxversion,inp,volfrac,rmin,penalty, volume_scaling, passivefilename):
    
    volConstraint=rwd.getVolconstr()
    total_volume = rwd.getTotalVol() # maximum volume of solid domain, vmax

    scaling = volume_scaling/(volfrac*total_volume) 

    if grad.size > 0:
        g1=rwd.getCurrentVolsens()

        if isinstance(passivefilename, str):
            passiveEl = rwd.readElemList(passivefilename)
            g1 = rwd.setElemVal(g1, passiveEl,0* np.ones(np.size(passiveEl))) 
        

        grad[:]=g1*scaling
    return volConstraint*scaling


def runCCX(ccx_version,inpFile,penalty,rmin,volfrac,itertop, nnz, q, beta, eta):
#run calxulix to get function and gradients

    a=subprocess.call([ccx_version, "-i", inpFile, "-p", str(penalty), "-r", str(rmin), "-v", str(volfrac), "-s", str(itertop),"-f", str(nnz), "-q", str(q), str(q), "-b", str(beta), "-e", str(eta)])