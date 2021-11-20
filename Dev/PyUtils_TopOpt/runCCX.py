#!/usr/bin/env python
from __future__ import print_function
import subprocess

def runCCX(ccx_version,inpFile,penalty,rmin,volfrac,itertop, nnz, q):
#run calxulix to get function and gradients

    a=subprocess.call([ccx_version, "-i", inpFile, "-p", str(penalty), "-r", str(rmin), "-v", str(volfrac), "-s", str(itertop),"-f", str(nnz), "-q", str(q)])
