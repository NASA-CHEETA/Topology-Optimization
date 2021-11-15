# Description :: This script imports the .msh file for node and element information.
# This script also calls the displacements.py module to import nodal discplacements in the solid
# Two write modes exist : 1) deform_mesh_str -- > Deforms a structured mesh
#                         2) deform_mesh_unstr ..> Deforms an unstructured mesh


import numpy as np
import pandas as pd
import math
import sys
import displacements

def deform_mesh_str(Mesh_name,Nnodes,Nelems):

    #Extract the node information
    Nodes = pd.read_csv(Mesh_name,header = None, sep = ',', skiprows = 2, nrows=Nnodes)

    #Extract the element information

    Element = pd.read_csv(Mesh_name,header = None, sep = ',', skiprows = Nnodes + 7, nrows=Nelems)

    LOC = Nodes.to_numpy()
    ELE = Element.to_numpy()

    Eleid = ELE[0:Nelems,0]
    E1 = ELE[0:Nelems,1]
    E2 = ELE[0:Nelems,2]
    E3 = ELE[0:Nelems,3]
    E4 = ELE[0:Nelems,4]
    E5 = ELE[0:Nelems,5]
    E6 = ELE[0:Nelems,6]
    E7 = ELE[0:Nelems,7]
    E8 = ELE[0:Nelems,8]

    Eleid= Eleid.reshape(Nelems,1)
    E1= E1.reshape(Nelems,1)
    E2= E2.reshape(Nelems,1)
    E3= E3.reshape(Nelems,1)
    E4= E4.reshape(Nelems,1)
    E5= E5.reshape(Nelems,1)
    E6= E6.reshape(Nelems,1)
    E7= E7.reshape(Nelems,1)
    E8= E8.reshape(Nelems,1)

    if Eleid[-1]!=Nelems:
        sys.exit("ERROR : Inconsistent number of elements extracted during mesh deformation")

    Nodeid = LOC[0:Nnodes,0]
    X = LOC[0:Nnodes,1]
    Y = LOC[0:Nnodes,2]
    Z = LOC[0:Nnodes,3]

    Nodeid= Nodeid.reshape(Nnodes,1)
    X = X.reshape(Nnodes,1)
    Y = Y.reshape(Nnodes,1)
    Z = Z.reshape(Nnodes,1)

    if Nodeid[-1]!=Nnodes:
        sys.exit("ERROR : Inconsistent number of nodes extracted during mesh deformation")

    U = displacements.displacements(Nnodes)

    NODEID_Cvrg = U[:,0]
    UX_Cvrg = U[:,1]
    UY_Cvrg = U[:,2]
    UZ_Cvrg = U[:,3]

    NODEID_Cvrg = NODEID_Cvrg.reshape(Nnodes,1)
    UX_Cvrg = UX_Cvrg.reshape(Nnodes,1)
    UY_Cvrg = UY_Cvrg.reshape(Nnodes,1)
    UZ_Cvrg = UZ_Cvrg.reshape(Nnodes,1)

    for i in range(Nnodes):
        if(NODEID_Cvrg[i,0]!=Nodeid[i,0]):
            sys.exit("ERROR : Displacement nodes not equal to mesh nodes")

    # Compute New Node location

    X = X + UX_Cvrg;
    Y = Y + UY_Cvrg;
    Z = Z + UZ_Cvrg;
    # Round the new nodal coordinates to 6 digits of precicion
    for i in range(Nnodes):
        X[i,0] = round(X[i,0],6)
        Y[i,0] = round(Y[i,0],6)
        Z[i,0] = round(Z[i,0],6)
    # Write New mesh file

    file = open("sample_deform.msh", "w")
    file.write("  **This containes mesh information")
    file.write("\n  *NODE, NSET=NALL")
    file.write("\n")
    for i in range(Nnodes):
        file.write("\n" '  ' + str(int(NODEID_Cvrg[i,0])) +'\t , \t'+ str(X[i,0]) +'  \t, \t  '+ str(Y[i,0]) + '  , ' + str(Z[i,0]))
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("*ELEMENT, TYPE=C3D8, ELSET=EALL")
    for i in range(Nelems):
        file.write("\n" + str(int(Eleid[i,0])) +' , '+ str(E1[i,0]) +' , '+ str(E2[i,0]) + ' , ' + str(E3[i,0])+ ' , ' +str(E4[i,0])+ ' , ' +str(E5[i,0])+ ' , '+ str(E6[i,0]) + ' , '+str(E7[i,0])+ ' , '+ str(E8[i,0]))
    file.write("\n")
    file.close()

    print(' Deformed STRUCTURED mesh succesfully written! ')


###########################################################################################
##########################################################################################
##########################################################################################


def deform_mesh_unstr(Mesh_name,Nnodes,Nelems):
    #Extract the node information
    Nodes = pd.read_csv(Mesh_name,header = None, sep = ',', skiprows = 2, nrows=Nnodes)

    #Extract the element information

    Element = pd.read_csv(Mesh_name,header = None, sep = ',', skiprows = Nnodes + 7, nrows=Nelems)

    LOC = Nodes.to_numpy()
    ELE = Element.to_numpy()

    Eleid = ELE[0:Nelems,0]
    E1 = ELE[0:Nelems,1]
    E2 = ELE[0:Nelems,2]
    E3 = ELE[0:Nelems,3]
    E4 = ELE[0:Nelems,4]


    Eleid= Eleid.reshape(Nelems,1)
    E1= E1.reshape(Nelems,1)
    E2= E2.reshape(Nelems,1)
    E3= E3.reshape(Nelems,1)
    E4= E4.reshape(Nelems,1)


    if Eleid[-1]!=Nelems:
        sys.exit("ERROR : Inconsistent number of elements extracted during mesh deformation")

    Nodeid = LOC[0:Nnodes,0]
    X = LOC[0:Nnodes,1]
    Y = LOC[0:Nnodes,2]
    Z = LOC[0:Nnodes,3]

    Nodeid= Nodeid.reshape(Nnodes,1)
    X = X.reshape(Nnodes,1)
    Y = Y.reshape(Nnodes,1)
    Z = Z.reshape(Nnodes,1)

    if Nodeid[-1]!=Nnodes:
        sys.exit("ERROR : Inconsistent number of nodes extracted during mesh deformation")

    U = displacements.displacements(Nnodes)

    NODEID_Cvrg = U[:,0]
    UX_Cvrg = U[:,1]
    UY_Cvrg = U[:,2]
    UZ_Cvrg = U[:,3]

    NODEID_Cvrg = NODEID_Cvrg.reshape(Nnodes,1)
    UX_Cvrg = UX_Cvrg.reshape(Nnodes,1)
    UY_Cvrg = UY_Cvrg.reshape(Nnodes,1)
    UZ_Cvrg = UZ_Cvrg.reshape(Nnodes,1)

    for i in range(Nnodes):
        if(NODEID_Cvrg[i,0]!=Nodeid[i,0]):
            sys.exit("ERROR : Displacement nodes not equal to mesh nodes")

    # Compute New Node location

    X = X + UX_Cvrg;
    Y = Y + UY_Cvrg;
    Z = Z + UZ_Cvrg;

    # Round the new nodal coordinates to 6 digits of precicion
    for i in range(Nnodes):
        X[i,0] = round(X[i,0],6)
        Y[i,0] = round(Y[i,0],6)
        Z[i,0] = round(Z[i,0],6)
    # Write New mesh file

    file = open("Solid/wing_deform.nam", "w")
    file.write("  **This containes mesh information")
    file.write("\n  *NODE, NSET=NALL")
    file.write("\n")
    for i in range(Nnodes):
        file.write("\n" '  ' + str(int(Nodeid[i,0])) +'\t , \t'+ str(X[i,0]) +'  \t, \t  '+ str(Y[i,0]) + '  , ' + str(Z[i,0]))
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("*ELEMENT, TYPE=C3D4, ELSET=EALL")
    for i in range(Nelems):
        file.write("\n" + str(int(Eleid[i,0])) +' , '+ str(E1[i,0]) +' , '+ str(E2[i,0]) + ' , ' + str(E3[i,0])+ ' , ' +str(E4[i,0]))
    file.write("\n")
    file.close()

    print(' Deformed UNSTRUCTURED mesh succesfully written! ')
