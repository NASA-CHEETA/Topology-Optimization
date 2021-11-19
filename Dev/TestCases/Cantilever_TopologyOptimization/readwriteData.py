#!/usr/bin/env python

from __future__ import print_function
#import os
import numpy as np


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
    Objfile="objectives.dat"
    Objs=np.loadtxt(Objfile,delimiter=',',dtype=float,usecols=0)
    v=float(Objs)
    return v


def getVolconstr():
    """
    Obtain total volume constraint value; rho*v-Volfrac*vinitial
    """
    Objfile="objectives.dat"

    Objs=np.loadtxt(Objfile,delimiter=',',dtype=float,usecols=1)
    v=float(Objs)
    return v

def getTotalVol():
    """
    Obtain total volume without density
    """
    Objfile="objectives.dat"

    Objs=np.loadtxt(Objfile,delimiter=',',dtype=float,usecols=2)
    return Objs

def getCurrentVol():
    """
    Obtain current volume based on filtered density
    """
    Objfile="objectives.dat"

    Objs=np.loadtxt(Objfile,delimiter=',',dtype=float,usecols=3)
    return Objs


def getComplsens():
    """
    Obtain filtered compliance sensitivity
    """
    Gradfile="sens_compliance.dat"

    Grads=np.loadtxt(Gradfile,delimiter=',',dtype=float,usecols=1)
    return Grads

def getCurrentVolsens():
    """
    Obtain filtered volume sensitivity
    """
    Gradfile="sens_volume.dat"

    Grads=np.loadtxt(Gradfile,delimiter=',',dtype=float,usecols=2)
    return Grads    

def getRhoPhys():
    """
    Obtain filtered density
    """
    file="rhos.dat"

    rho=np.loadtxt(file,delimiter=',',dtype=float,usecols=1)
    return rho


def writeCurrentrho(rho):
    """
    write current unfiltered rho to file
    """
    rhofile="density.dat"

    np.savetxt(rhofile,rho,fmt='%1.6f')

def writerhoPhys(rhoPhys):
    """
    write current filtered rho to file
    """
    rhofile="rhoPhys.dat"

    np.savetxt(rhofile,rhoPhys,fmt='%1.6f')


def postprocess(coord,map,outputinp,rhoPhys,cutoff):
    outputfile="mapelementsout.msh"

    with open(map, "r") as f:
            lines = f.readlines()
    with open(outputfile, "w") as new_f:
        new_f.write(lines[0])
        for line in lines[1:]:  
            b=line.split(",")
               # print(int(b[0]))
            if rhoPhys[int(b[0])-1] > cutoff :       
                new_f.write(line)

    l1="*HEADING\n"
    l2="Optimized output\n\n"

    l3="*INCLUDE,INPUT="+coord
    l4="\n*INCLUDE,INPUT="+outputfile

    with open(outputinp, "w") as text_file:
            text_file.writelines([l1,l2,l3,l4])
            print("writing final output done!")