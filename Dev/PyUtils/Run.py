# Description : This is the main execution file for launching aeroelastic optimziation.
# This script calls the shell script for the partitioned aero-elastic solver
# Upon convergence and an Exit Code 0, the deform mesh module is called to create a deformed mesh.at aero-elastic equillibrium
import numpy as np
import pandas as pd
import math
import FSI
import deform_mesh as dm
import sys

# Node and element count is requred for parsing the forces and .msh file for
# mesh deformation.
Nnodes = 10467
Nelems = 45513
Mesh_type = 0   #Unstructured mesh
Mesh_name = 'Solid/wing3dmsh.nam'

run_FSI = False

if run_FSI == True:
    Flag = FSI.AeroElastic()

    if ( Flag != 0):
        sys.exit("HALT!: Aeroelastic Simulations Diverged! ")
else:
    print("Running in debug mode! ")

if Mesh_type ==0:
    dm.deform_mesh_unstr(Mesh_name,Nnodes,Nelems)

else:
    dm.deform_mesh_str(Mesh_name,Nnodes,Nelems)
