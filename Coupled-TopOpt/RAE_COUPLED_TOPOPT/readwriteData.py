#!/usr/bin/env python

from __future__ import print_function
#import os
import numpy as np


def writeHeader(ccxversion):
    print("\t\t\t*******************************************************************")
    print("\t\t\t\t\tTOPOLOGY OPTIMIZATION USING " + ccxversion+"\t\t")
    print("\t\t\t\t\tCDILAB, Dept of Aerospace Engineering, UIUC\t\t")
    print("\t\t\t*******************************************************************")


def maxnnz():
    """
    Find the maximum non zero in the density matrix from the dnnz.dat file
    """
    nnz = np.loadtxt("dnnz.dat")
    return np.max(nnz)


def getCompl():
    """
    Obtain compliance value
    """
    Objfile = "objectives.dat"
    Objs = np.loadtxt(Objfile, delimiter=',', dtype=float, usecols=0)
    v = float(Objs)
    return v


def getVolconstr():
    """
    Obtain total volume constraint value; rho*v-Volfrac*vinitial
    """
    Objfile = "objectives.dat"
    Objs = np.loadtxt(Objfile, delimiter=',', dtype=float, usecols=1)
    v = float(Objs)
    return v


def getTotalVol():
    """
    Obtain total volume without density
    """
    Objfile = "objectives.dat"
    Objs = np.loadtxt(Objfile, delimiter=',', dtype=float, usecols=2)
    return Objs


def getCurrentVol():
    """
    Obtain current volume based on filtered density
    """
    Objfile = "objectives.dat"
    Objs = np.loadtxt(Objfile, delimiter=',', dtype=float, usecols=3)
    return Objs


def rho_column(file, column=0):
    """
    Obtain a column from a comma separated file
    """
    rho = np.loadtxt(file, delimiter=',', dtype=float, usecols=column)
    return rho


def getComplsens():
    """
    Obtain filtered compliance sensitivity
    """
    Gradfile = "sens_compliance.dat"
    Grads = np.loadtxt(Gradfile, delimiter=',', dtype=float, usecols=1)
    return Grads


def getCurrentVolsens():
    """
    Obtain filtered volume sensitivity
    """
    Gradfile = "sens_volume.dat"
    Grads = np.loadtxt(Gradfile, delimiter=',', dtype=float, usecols=2)
    return Grads


def getRhoPhys():
    """
    Obtain filtered density
    """
    file = "rhos.dat"
    rho = np.loadtxt(file, delimiter=',', dtype=float, usecols=1)
    return rho


def writeCurrentrho(rho):
    """
    write current unfiltered rho to file
    """
    rhofile = "density.dat"
    np.savetxt(rhofile, rho, fmt='%1.6f')


def writerhoPhys(rhoPhys):
    """
    write current filtered rho to file
    """
    rhofile = "density_FSI.dat"
    np.savetxt(rhofile, rhoPhys, fmt='%1.6f')


def postprocess(coord, map, outputinp, rhoPhys, cutoff):
    outputfile = "mapelementsout.msh"

    with open(map, "r") as f:
        lines = f.readlines()
    with open(outputfile, "w") as new_f:
        new_f.write(lines[0])
        for line in lines[1:]:
            b = line.split(",")
            # print(int(b[0]))
            if rhoPhys[int(b[0])-1] > cutoff:
                new_f.write(line)

    l1 = "*HEADING\n"
    l2 = "Optimized output\n\n"

    l3 = "*INCLUDE,INPUT="+coord
    l4 = "\n*INCLUDE,INPUT="+outputfile

    with open(outputinp, "w") as text_file:
        text_file.writelines([l1, l2, l3, l4])
        print("writing final output done!")


def parse_fea_set(File, outputfile):
    """
    Reads the ccx set file, and writes into a file as 
    a column vector without comma
    """
    # Read the file
    with open(File, "r") as f:
        lines = f.readlines()
        # Write to file
        with open(outputfile, "w") as new_f:
            for line in lines[1:]:
                b = line.split(",")
                for each in b[:-1]:
                    new_f.write(each)
                    new_f.write('\n')
            # The last term in last line was not comma, but ignored.
            # so write it
            new_f.write(b[-1])


def readElemList(File='passiveElemList', col=0):
    """
    Read a column of element numbers from a file. Use this to set void/passive elements
    """
    el = np.loadtxt(File, dtype=int)
    return el


def setElemVal(vector, index, value):
    """
    Set entries at index locations in vector to a value.
    Use it to assign 0 or 1 values. 
    """
    vector[index-1] = value
    return vector


def exclude_elements(excludefilenameList, x0, value):
    """
    Sets all entries in vector x0 to value for all files in list ['f1.nam', 'f2,nam', ...]
    """
    for eachFile in excludefilenameList:
        if isinstance(eachFile, str):
            #print('Modifying entries from ' + eachFile)
            passiveEl = readElemList(eachFile)
            x0 = setElemVal(x0, passiveEl, value * np.ones(np.size(passiveEl)))
    return x0


def update_beta(iter, iter_init, iter_beta, beta, beta_factor, betamax):
    """
    Update new value of beta, Based on https://github.com/thsmit/TopOpt_in_PETSc_wrapped_in_Python/blob/master/topoptlib/src/Filter.cc
    """
    #if iter>15 and iter<=iter_init and beta<=betamax:
    #   beta = beta+1

    #el
    if iter >= iter_init and iter % iter_beta == 0 and beta <= betamax:
        beta = beta*beta_factor
    if beta > betamax:
        beta = betamax
    return beta

    ##############################################################################


def apply_heavisidefilter(x, eta, beta):
    """
    Apply tanh filter to rho, also calculate derivative
    """
    f = 0 + (np.tanh(beta*eta) + np.tanh(beta*(x-eta))) / \
        (np.tanh(beta*eta) + np.tanh(beta*(1-eta)))
    df = beta*(1-(np.tanh(beta*(x-eta)))**2) / \
        (np.tanh(beta*eta) + np.tanh(beta*(1-eta)))
    return f, df


def sum_projection(x, eta, beta):
    """
    Calculate difference in volume by heaviside operator
    """
    x_sum = np.sum(x)
    xh = apply_heavisidefilter(x, eta, beta)[0]
    xh_sum = np.sum(xh)
    return x_sum - xh_sum


def my_bisection(f, a, b, tol, x, beta):
    # approximates a root, R, of f bounded
    # by a and b to within tolerance
    # | f(m) | < tol with m the midpoint
    # between a and b Recursive implementation
    # get midpoint
    m = (a + b)/2
    if np.abs(f(x, m, beta)) < tol:
        # stopping condition, report m as root
        return m
    elif f(x, a, beta) * f(x, m, beta) > 0:
        # case where m is an improvement on a.
        # Make recursive call with a = m
        return my_bisection(f, m, b, tol, x, beta)
    elif f(x, b, beta)*f(x, m, beta) > 0:
        # case where m is an improvement on b.
        # Make recursive call with b = m
        return my_bisection(f, a, m, tol, x, beta)


def update_eta(x, beta):
    """
    Determine volume preserving eta using bisection search
    """
    a = 0
    b = 1
    tol = 1e-4
    return my_bisection(sum_projection, a, b, tol, x, beta)
