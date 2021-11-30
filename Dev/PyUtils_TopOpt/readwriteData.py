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

def rho_column(file, column = 0):
    """
    Obtain a column from a comma separated file
    """
    rho=np.loadtxt(file,delimiter=',',dtype=float,usecols=column)
    return rho   


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
    rhofile="density_FSI.dat"

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
                b=line.split(",")    
                for each in b[:-1]:
                    new_f.write(each)
                    new_f.write('\n')
            # The last term in last line was not comma, but ignored.
            # so write it
            new_f.write(b[-1])


def read_fea_set(File):
    """
    Read a column of element numbers or files written by 
    parse_fea_set(). Use it to set active/passive elements
    """
    el=np.loadtxt(File, delimiter='\n',dtype=int)
    return el

def value_fea_set(vector, index, value):
    """
    Set entries at index locations in vector to a value.
    Use it to assign 0 or 1 values. 
    """
    vector[index-1] = value
    return vector