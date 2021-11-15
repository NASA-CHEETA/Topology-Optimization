import numpy as np
import pandas as pd
import sys

def interp_tet(mesh,inputrhofile, Nnodes, Nelems):
    """
    To interpolate tetrahedral element density to nodes from rhos.dat obtained from ccx.
    Requires node and element information.

    #Ghanendra K Das, UIUC, 2021
    """
    #Extract the node information
    Nodes = pd.read_csv(mesh,header = None, sep = ',', skiprows = 2, nrows=Nnodes)

    #Extract the element information
    Element = pd.read_csv(mesh,header = None, sep = ',', skiprows = Nnodes+2+2, nrows=Nelems)

    #Extract the density information
    Rhos = pd.read_csv(inputrhofile,header = None, sep = ',', nrows=Nelems)

    LOC = Nodes.to_numpy()
    ELE = Element.to_numpy()
    RHO = Rhos.to_numpy()

    RHO_FILTERED=RHO[0:Nelems,1]        
    RHO_FILTERED=RHO_FILTERED.reshape(Nelems,1)
    RHO_UNFILTERED=RHO[0:Nelems,0]
    RHO_UNFILTERED=RHO_UNFILTERED.reshape(Nelems,1)
    densityinterp = np.zeros(Nnodes)    #To store sum of connecting element rhos at node

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
        sys.exit("ERROR : Inconsistent number of elements read from mesh file")

    Nodeid = LOC[0:Nnodes,0]
    X = LOC[0:Nnodes,1]
    Y = LOC[0:Nnodes,2]
    Z = LOC[0:Nnodes,3]

    Nodeid= Nodeid.reshape(Nnodes,1)
    X = X.reshape(Nnodes,1)
    Y = Y.reshape(Nnodes,1)
    Z = Z.reshape(Nnodes,1)

    if Nodeid[-1]!=Nnodes:
        sys.exit("ERROR : Inconsistent number of nodes read from mesh file")

    #-----------------------------------------------------------------------------#
    ## Extract nodal density
    #-----------------------------------------------------------------------------#
    print("Finding node rho from elements...")
    #Look for all elements connecting at a node, and add their rhos
    for elcount, elem in enumerate(Eleid):
        rhoe=RHO_FILTERED[elem-1]

        densityinterp[E1[elem-1]-1]= max(densityinterp[E1[elem-1]-1],rhoe)
        densityinterp[E2[elem-1]-1]= max(densityinterp[E2[elem-1]-1],rhoe)
        densityinterp[E3[elem-1]-1]= max(densityinterp[E3[elem-1]-1],rhoe)
        densityinterp[E4[elem-1]-1]= max(densityinterp[E4[elem-1]-1],rhoe)

    #Average density at node
    densityinterp_ave = densityinterp.reshape(Nnodes,1)

    return Nodeid, X, Y, Z, Eleid, E1, E2, E3, E4, densityinterp_ave 
    

if __name__ == "__main__":
    mesh ="wing_deform.msh"
    inputrhofile = "rhos.dat"
 
    Nnodes= 10467
    Nelems= 45513

    outputfile= "rhoInterpolated.dat"


    Nodeid, X, Y, Z, Eleid, E1, E2, E3, E4, densityinterp_ave =interp_tet(mesh,inputrhofile, Nnodes, Nelems)

    file = open(outputfile, "w")
    file.write("TITLE = \"Visualization of the solid solution\"")
    file.write("\n")
    file.write("VARIABLES = \"x\"\"y\"\"z\"\"Rho\"")
    file.write("\n")
    file.write("ZONE STRANDID=2, SOLUTIONTIME=1, NODES="+str(Nnodes)+", ELEMENTS="+ str(Nelems)+", DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON")
    for i in range(Nnodes):
        file.write("\n" +"{:.5e}" .format(X[i,0]) +'  \t  '+ "{:.5e}".format(Y[i,0]) + ' \t \t ' + "{:.5e}".format(Z[i,0]) + ' \t ' +"{:.5e}".format(densityinterp_ave[i,0]))

    for i in range(Nelems):
        file.write("\n" + str(E1[i,0]) +' \t '+ str(E2[i,0]) + ' \t ' + str(E3[i,0])+ ' \t ' +str(E4[i,0]))
    file.write("\n")
    file.close()
    print("Complete!")