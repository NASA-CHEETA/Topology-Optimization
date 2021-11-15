import numpy as np
import sys
import os

def displacements(Nnodes):
    C = []
    print('Reading nodal displacements for mesh deformation ...')
    f = open('Solid/wing3d.csv','r')
    c = f.readlines()
    for i in range(len(c)):
        text = c[i].split( )
        try:
            C.append ([float(item) for item in text])
        except:
            continue
    C[:] = [x for x in C if x != []]
    nodeid = []
    ux = []
    uy = []
    uz = []
    for i in range(len(C)):
        nodeid.append(C[i][0])
        ux.append(C[i][1])
        uy.append(C[i][2])
        uz.append(C[i][3])
# Convert lists to Numpy Arrays

    #print(len(ux))
    #print(len(uy))
    #print(len(uz))

    NODEID = np.array(nodeid)
    UX = np.array(ux)
    UY = np.array(uy)
    UZ = np.array(uz)

# Reshape the arrays to form two-dimensional arrays
    NODEID = NODEID.reshape(NODEID.size,1)
    UX = UX.reshape(NODEID.size,1)
    UY = UY.reshape(NODEID.size,1)
    UZ = UZ.reshape(NODEID.size,1)

# Check if correct number of nodes are extracted

    if (NODEID[-1,0]!=Nnodes):
        sys.exit("ERROR : Inconsistent number of nodes during displacement analysis")
# Extract displacements from the last time step only
    NODEID_Cvrg = NODEID[:Nnodes,0]
    UX_Cvrg = UX[-Nnodes:,0]
    UY_Cvrg = UY[-Nnodes:,0]
    UZ_Cvrg = UZ[-Nnodes:,0]
    NODEID_Cvrg = NODEID_Cvrg.reshape(Nnodes,1)
    UX_Cvrg = UX_Cvrg.reshape(Nnodes,1)
    UY_Cvrg = UY_Cvrg.reshape(Nnodes,1)
    UZ_Cvrg = UZ_Cvrg.reshape(Nnodes,1)


    U = np.concatenate((NODEID_Cvrg,UX_Cvrg,UY_Cvrg,UZ_Cvrg),axis = 1)
    return U
