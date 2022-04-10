#!/usr/bin/env python

from __future__ import print_function
#import os
import numpy as np
import subprocess
import readwriteData as rwd


def callrunccx(x, grad, ccxversion, inp, volfrac, rmin, penalty, compfile,  plotIterationWiseFlag, scaling, passivefilename, voidfilename):
    """
    Call ccx and calculate compliance
    """

    # Update iteration counter
    callrunccx.counter += 1
    itertop = callrunccx.counter

    # Update filter order
    if (callrunccx.counter) % (callrunccx.q_iter) == 0 and ((callrunccx.q) < (callrunccx.q_max)):
        callrunccx.q *= callrunccx.q_factor

    # Update beta for Heaviside
    callrunccx.beta = rwd.update_beta(callrunccx.counter, callrunccx.beta_iter_init,
                                      callrunccx.beta_iter, callrunccx.beta, callrunccx.beta_factor, callrunccx.betamax)

    # Update dnnz values for i=2 onwards
    if callrunccx.counter == 2:
        temp = rwd.maxnnz()
        callrunccx.nnz = int(temp)

    # Write design variables to density.dat
    rwd.writeCurrentrho(x)

    # Call ccx to write densityfsi.dat
    """ 
    gkdas2: Prateek, Following call to ccx with caltop_mode=2 will filter the design variable, write filtered design
                    variables to densityfsi.dat and exit. Copy this file into the folder where you will perform
                    fsi analysis. This helps to avoid having to hold a giant density filter matrix in memory while performing
                    FSI analysis.

                    This call in very first iteration will also write drow.dat, dcol.dat, dval.dat and dnnz.dat.
                    They contain compressed density filter matrix and will be reused later in all subsequent caltop calls,
                    assuming that the cost of moving files will be less than actually recalculating them. 
    """
    caltop_mode = 2  # For filter mode. DO NOT CHANGE
    runCCX(ccxversion, inp, penalty, rmin, volfrac, itertop,
           callrunccx.nnz, callrunccx.q, caltop_mode, callrunccx.beta, callrunccx.eta)

    # You have the densityfsi.dat. Now perform the fsi with this densities. Finger crossed!

    # FSI has converged by now. Use the values to get compliance and sensitivities. Call ccx in full mode
    caltop_mode = 1  # For standard mode. DO NOT CHANGE
    runCCX(ccxversion, inp, penalty, rmin, volfrac, itertop,
           callrunccx.nnz, callrunccx.q, caltop_mode, callrunccx.beta, callrunccx.eta)

    # read rhoPhys and store iterationwise
    if plotIterationWiseFlag:
        rhoPhys = rwd.getRhoPhys()
        iterwiseRhofileName = "rho"+str(itertop)+".dat"
        np.savetxt(iterwiseRhofileName, rhoPhys, fmt='%1.6f')

    # Get compliance value
    compliance = rwd.getCompl()
    print('Compliance=')
    print(compliance)

    f = open(compfile, "a")
    f.write("\n")
    f.write(str(compliance))
    f.close()

    if grad.size > 0:
        g = rwd.getComplsens()

        # Supress sensitivity to passive and void elements
        g = rwd.exclude_elements(passivefilename, g, 0)
        g = rwd.exclude_elements(voidfilename, g, 0)

        grad[:] = g*scaling
    return compliance*scaling


def getVolconstraint(x, grad, ccxversion, inp, volfrac, rmin, penalty, volume_scaling, passivefilename, voidfilename):

    volConstraint = rwd.getVolconstr()
    total_volume = rwd.getTotalVol()  # maximum volume of solid domain, vmax

    #scaling = volume_scaling/(volfrac*total_volume)
    scaling = volume_scaling
    if grad.size > 0:
        g1 = rwd.getCurrentVolsens()

        g1 = rwd.exclude_elements(passivefilename, g1, 0)
        g1 = rwd.exclude_elements(voidfilename, g1, 0)

        grad[:] = g1*scaling
    return volConstraint*scaling


def runCCX(ccx_version, inpFile, penalty, rmin, volfrac, itertop, nnz, q, caltop_mode, beta, eta):
    # run calxulix to get function and gradients

    a = subprocess.call([ccx_version, "-i", inpFile, "-p", str(penalty), "-r", str(rmin), "-v", str(volfrac), "-s",
                        str(itertop), "-f", str(nnz), "-q", str(q), str(q), "-m", str(caltop_mode), "-b", str(beta), "-e", str(eta)])
