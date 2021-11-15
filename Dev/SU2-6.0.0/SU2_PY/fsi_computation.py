#!/usr/bin/env python

## \file fsi_computation.py
#  \brief Python wrapper code for FSI computation by coupling a third-party structural solver to SU2.
#  \author David Thomas
#  \version 6.0.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import time as timer
from math import *	# use mathematical expressions
from optparse import OptionParser	# use a parser for configuration

import SU2	# imports SU2 python tools
import FSI	# imports FSI python tools


# imports the CFD (SU2) module for FSI computation
import pysu2

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

  # --- Get the FSI conig file name form the command line options --- #
  parser=OptionParser()
  parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
  parser.add_option("--parallel", action="store_true",
                      help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

  (options, args)=parser.parse_args()

  if options.with_MPI == True:
    from mpi4py import MPI  # MPI is initialized from now by python and can be continued in C++ !
    comm = MPI.COMM_WORLD
    myid = comm.Get_rank()
    numberPart = comm.Get_size()
    have_MPI = True
  else:
    comm = 0
    myid = 0
    numberPart = 1
    have_MPI = False

  rootProcess = 0

  # --- Set the working directory --- #
  if myid == rootProcess:
      if os.getcwd() not in sys.path:
          sys.path.append(os.getcwd())
	  print("Setting working directory : {}".format(os.getcwd()))
      else: 
	  print("Working directory is set to {}".format(os.getcwd()))

  # starts timer
  start = timer.time()

  confFile = str(options.filename)

  FSI_config = FSI.io.FSIConfig(confFile) 		# FSI configuration file
  CFD_ConFile = FSI_config['CFD_CONFIG_FILE_NAME']	# CFD configuration file
  CSD_ConFile = FSI_config['CSD_CONFIG_FILE_NAME']	# CSD configuration file

  CSD_Solver = FSI_config['CSD_SOLVER']			# CSD solver

  if have_MPI == True:
    comm.barrier()

  # --- Initialize the fluid solver --- #
  if myid == rootProcess:
    print('\n***************************** Initializing fluid solver *****************************')
  try:
    FluidSolver = pysu2.CFluidDriver(CFD_ConFile, 1, FSI_config['NDIM'], comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CSingleZoneDriver : ',exception)
    if have_MPI == True:
      print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
    else:
      print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    return

  if have_MPI == True:
    comm.barrier()
  
  # --- Initialize the solid solver --- # (!! for now we are using only serial solid solvers)
  if myid == rootProcess:
    print('\n***************************** Initializing solid solver *****************************')
    if CSD_Solver == 'METAFOR':
      from MetaforSolver import MtfSolver
      SolidSolver = MtfSolver(CSD_ConFile)
    elif CSD_Solver == 'NATIVE':
      import NativeSolid
      SolidSolver = NativeSolid.NativeSolidSolver(CSD_ConFile, True)
    elif CSD_Solver == 'GETDP':
      import GetDPSolver
      SolidSolver = GetDPSolver.GetDPSolver(CSD_ConFile, True)
    elif CSD_Solver == 'TESTER':
      SolidSolver = FSI.PitchPlungeAirfoilStructuralTester.Solver(CSD_ConFile)
  else:
    SolidSolver = None

  if have_MPI == True:
    comm.barrier()

  # --- Initialize and set the FSI interface (coupling environement) --- #
  if myid == rootProcess:
    print('\n***************************** Initializing FSI interface *****************************')
  if have_MPI == True:
    comm.barrier()
  FSIInterface = FSI.Interface(FSI_config, FluidSolver, SolidSolver, have_MPI)
  
  if myid == rootProcess:
    print('\n***************************** Connect fluid and solid solvers *****************************')
  if have_MPI == True:
    comm.barrier()
  FSIInterface.connect(FSI_config, FluidSolver, SolidSolver)

  if myid == rootProcess:
    print('\n***************************** Mapping fluid-solid interfaces *****************************')
  if have_MPI == True:
    comm.barrier()
  FSIInterface.interfaceMapping(FluidSolver, SolidSolver, FSI_config)
 
  if have_MPI == True: 
    comm.barrier()

  # --- Launch a steady or unsteady FSI computation --- #
  if FSI_config['UNSTEADY_SIMULATION'] == "YES":
    try:
      FSIInterface.UnsteadyFSI(FSI_config, FluidSolver, SolidSolver)
    except NameError as exception:
      if myid == rootProcess:
        print('An NameError occured in FSIInterface.UnsteadyFSI : ',exception)
    except TypeError as exception:
      if myid == rootProcess:
        print('A TypeError occured in FSIInterface.UnsteadyFSI : ',exception)
    except KeyboardInterrupt as exception :
      if myid == rootProcess:
        print('A KeyboardInterrupt occured in FSIInterface.UnsteadyFSI : ',exception)
  else:
    try:
      NbExtIter = FSI_config['NB_EXT_ITER']
      FSIInterface.SteadyFSI(FSI_config, FluidSolver, SolidSolver)
    except NameError as exception:
      if myid == rootProcess:
        print('An NameError occured in FSIInterface.SteadyFSI : ',exception)
    except TypeError as exception:
      if myid == rootProcess:
        print('A TypeError occured in FSIInterface.SteadyFSI : ',exception)
    except KeyboardInterrupt as exception :
      if myid == rootProcess:
        print('A KeyboardInterrupt occured in FSIInterface.SteadyFSI : ',exception)
  
  if have_MPI == True:
    comm.barrier()

  # --- Exit cleanly the fluid and solid solvers --- #
  FluidSolver.Postprocessing()
  if myid == rootProcess:
      SolidSolver.exit()

  if have_MPI == True:
    comm.barrier()

  # stops timer
  stop = timer.time()
  elapsedTime = stop-start
  
  if myid == rootProcess:
    print("\n Computation successfully performed in {} seconds.".format(elapsedTime))

  return

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
