import numpy as np
import pandas as pd
import sys

def interp_hybrid(mesh,rho,Nnodes,Nelems, Nelemshex,Nelemstet):
    """
    To interpolate element density to nodes from rhos.dat obtained from ccx.
    Requires node and element information. This version has tet+hex capability
    The value at a node = max rho among all connecting elements at the node
    #Ghanendra K Das, UIUC, 2021

    Inputs: mesh file, element rho, total nodes, total elements, number of hex elements, number of tet elements.
    Note: Hex elements must be listed ahead of tet element in the mesh file

    TO DO: multiple element sets
    """
    #Extract the node information
    Nodes = pd.read_csv(mesh,header = None, sep = ',', skiprows = 2, nrows=Nnodes)

    #Extract the element information for hex
    Elementhex = pd.read_csv(mesh,header = None, sep = ',', skiprows = Nnodes+2+2, nrows=Nelemshex)

    #Extract the element information for hex
    Elementtet = pd.read_csv(mesh,header = None, sep = ',', skiprows = Nnodes+Nelemshex+2+2+2, nrows=Nelemstet)

    #Extract the density information
    Rhos = pd.read_csv(rho,header = None, sep = ',', nrows=Nelems)

    LOC = Nodes.to_numpy()
    ELEHEX = Elementhex.to_numpy()
    ELETET = Elementtet.to_numpy()    
    RHO = Rhos.to_numpy()

    RHO_FILTERED=RHO[0:Nelems,1]        
    RHO_FILTERED=RHO_FILTERED.reshape(Nelems,1)
    RHO_UNFILTERED=RHO[0:Nelems,0]
    RHO_UNFILTERED=RHO_UNFILTERED.reshape(Nelems,1)
    densityinterp = np.zeros(Nnodes)    #To store sum of connecting element rhos at node
    el_to_node_count = np.zeros(Nnodes) #To store number of elements sharing each node

    EleidHEX = ELEHEX[0:Nelemshex,0]
    E1HEX =ELEHEX[0:Nelemshex,1]
    E2HEX = ELEHEX[0:Nelemshex,2]
    E3HEX = ELEHEX[0:Nelemshex,3]
    E4HEX = ELEHEX[0:Nelemshex,4]
    E5HEX = ELEHEX[0:Nelemshex,5]
    E6HEX = ELEHEX[0:Nelemshex,6]
    E7HEX = ELEHEX[0:Nelemshex,7]
    E8HEX = ELEHEX[0:Nelemshex,8]

    EleidHEX= EleidHEX.reshape(Nelemshex,1)
    E1HEX= E1HEX.reshape(Nelemshex,1)
    E2HEX= E2HEX.reshape(Nelemshex,1)
    E3HEX= E3HEX.reshape(Nelemshex,1)
    E4HEX= E4HEX.reshape(Nelemshex,1)
    E5HEX= E5HEX.reshape(Nelemshex,1)
    E6HEX= E6HEX.reshape(Nelemshex,1)
    E7HEX= E7HEX.reshape(Nelemshex,1)
    E8HEX= E8HEX.reshape(Nelemshex,1)


    EleidTET = ELETET[0:Nelemstet,0]
    E1TET = ELETET[0:Nelemstet,1]
    E2TET = ELETET[0:Nelemstet,2]
    E3TET = ELETET[0:Nelemstet,3]
    E4TET = ELETET[0:Nelemstet,4]

    EleidTET= EleidTET.reshape(Nelemstet,1)
    E1TET= E1TET.reshape(Nelemstet,1)
    E2TET= E2TET.reshape(Nelemstet,1)
    E3TET= E3TET.reshape(Nelemstet,1)
    E4TET= E4TET.reshape(Nelemstet,1)


    if Nelems!=(Nelemstet+Nelemshex):
        sys.exit("ERROR : Inconsistent number of elements read from mesh file. Number of hex and tets do not add up to total number of elements. \n")   

    Nodeid = LOC[0:Nnodes,0]
    X = LOC[0:Nnodes,1]
    Y = LOC[0:Nnodes,2]
    Z = LOC[0:Nnodes,3]

    Nodeid= Nodeid.reshape(Nnodes,1)
    X = X.reshape(Nnodes,1)
    Y = Y.reshape(Nnodes,1)
    Z = Z.reshape(Nnodes,1)


    #-----------------------------------------------------------------------------#
    ## Extract nodal density
    #-----------------------------------------------------------------------------#
    print("Determining values wrt tet elements...")
    #Look for all elements connecting at a node, and add their rhos in tets
    for elcount, elem in enumerate(EleidTET):
        rhoe=RHO_FILTERED[elem-1]
    
        elemtet=elem-Nelemshex

        densityinterp[E1TET[elemtet-1]-1]= max(densityinterp[E1TET[elemtet-1]-1],rhoe)
        densityinterp[E2TET[elemtet-1]-1]= max(densityinterp[E2TET[elemtet-1]-1],rhoe)
        densityinterp[E3TET[elemtet-1]-1]= max(densityinterp[E3TET[elemtet-1]-1],rhoe)
        densityinterp[E4TET[elemtet-1]-1]= max(densityinterp[E4TET[elemtet-1]-1],rhoe)

    print("Determining values wrt hex elements...")
    #Look for all elements connecting at a node, and add their rhos in hex
    for elcount, elem in enumerate(EleidHEX):
     
        rhoe=RHO_FILTERED[elem-1]
        
        densityinterp[E1HEX[elem-1]-1]= max(densityinterp[E1HEX[elem-1]-1],rhoe)
        densityinterp[E2HEX[elem-1]-1]= max(densityinterp[E2HEX[elem-1]-1],rhoe)
        densityinterp[E3HEX[elem-1]-1]= max(densityinterp[E3HEX[elem-1]-1],rhoe)
        densityinterp[E4HEX[elem-1]-1]= max(densityinterp[E4HEX[elem-1]-1],rhoe)
        densityinterp[E5HEX[elem-1]-1]= max(densityinterp[E5HEX[elem-1]-1],rhoe)
        densityinterp[E6HEX[elem-1]-1]= max(densityinterp[E6HEX[elem-1]-1],rhoe)
        densityinterp[E7HEX[elem-1]-1]= max(densityinterp[E7HEX[elem-1]-1],rhoe)
        densityinterp[E8HEX[elem-1]-1]= max(densityinterp[E8HEX[elem-1]-1],rhoe)

    #Density at nodes density at node
    densityinterp_ave = densityinterp.reshape(Nnodes,1)
    return Nodeid, X,Y,Z, EleidHEX, E1HEX, E2HEX, E3HEX, E4HEX, E5HEX, E6HEX, E7HEX, E8HEX, EleidTET, E1TET, E2TET, E3TET, E4TET, densityinterp_ave

if __name__ == "__main__":
    mesh = "wing3Dmsh.nam"
    rho = "rhosq15.dat"
    outputfile= "combined.dat"
    Nnodes= 46882
    Nelems= 99368
    Nelemstet = 70849
    Nelemshex =  28519
    Nodeid, X,Y,Z, EleidHEX, E1HEX, E2HEX, E3HEX, E4HEX, E5HEX, E6HEX, E7HEX, E8HEX, EleidTET, E1TET, E2TET, E3TET, E4TET, densityinterp_ave = interp_hybrid(mesh,rho,Nnodes,Nelems, Nelemshex,Nelemstet)

    print('Writing to a file')

    '''
    file = open("hexalone.dat", "w")
    file.write("TITLE = \"Visualization of the solid solution\"")
    file.write("\n")
    file.write("VARIABLES = \"x\"\"y\"\"z\"\"Rho\"")
    file.write("\n")
    file.write("ZONE STRANDID=2, SOLUTIONTIME=1, NODES="+str(Nnodes)+", ELEMENTS="+ str(Nelemshex)+", DATAPACKING=POINT, ZONETYPE=FEBRICK")
    for i in range(Nnodes):
        file.write("\n" +"{:.5e}" .format(X[i,0]) +'  \t  '+ "{:.5e}".format(Y[i,0]) + ' \t \t ' + "{:.5e}".format(Z[i,0]) + ' \t ' +"{:.5e}".format(densityinterp_ave[i,0]))

    for i in range(Nelemshex):
        file.write("\n" + str(E1HEX[i,0]) +' \t '+ str(E2HEX[i,0]) + ' \t ' + str(E3HEX[i,0])+ ' \t ' +str(E4HEX[i,0]) +' \t '+str(E5HEX[i,0]) +' \t '+ str(E6HEX[i,0]) + ' \t ' + str(E7HEX[i,0])+ ' \t ' +str(E8HEX[i,0]))
    file.write("\n")
    file.close()


    #######################################
    file = open("tetalone.dat", "w")

    file.write("TITLE = \"Visualization of the solid solution\"")
    file.write("\n")
    file.write("VARIABLES = \"x\"\"y\"\"z\"\"Rho\"")
    file.write("\n")
    file.write("ZONE STRANDID=2, SOLUTIONTIME=1, NODES="+str(Nnodes)+", ELEMENTS="+ str(Nelemstet)+", DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON")
    for i in range(Nnodes):
        file.write("\n" +"{:.5e}" .format(X[i,0]) +'  \t  '+ "{:.5e}".format(Y[i,0]) + ' \t \t ' + "{:.5e}".format(Z[i,0]) + ' \t ' +"{:.5e}".format(densityinterp_ave[i,0]))

    for i in range(Nelemstet):
        file.write("\n" + str(E1TET[i,0]) +' \t '+ str(E2TET[i,0]) + ' \t ' + str(E3TET[i,0])+ ' \t ' +str(E4TET[i,0]))
    file.write("\n")
    file.close()

    #####################################################################
    '''


    file = open(outputfile, "w")

    file.write("TITLE = \"Visualization of the solid solution\"")
    file.write("\n")
    file.write("VARIABLES = \"x\"\"y\"\"z\"\"Rho\"")
    file.write("\n")
    file.write("ZONE STRANDID=2, SOLUTIONTIME=1, NODES="+str(Nnodes)+", ELEMENTS="+ str(Nelemshex)+", DATAPACKING=POINT, ZONETYPE=FEBRICK")
    for i in range(Nnodes):
        file.write("\n" +"{:.5e}" .format(X[i,0]) +'  \t  '+ "{:.5e}".format(Y[i,0]) + ' \t \t ' + "{:.5e}".format(Z[i,0]) + ' \t ' +"{:.5e}".format(densityinterp_ave[i,0]))

    for i in range(Nelemshex):
        file.write("\n" + str(E1HEX[i,0]) +' \t '+ str(E2HEX[i,0]) + ' \t ' + str(E3HEX[i,0])+ ' \t ' +str(E4HEX[i,0]) +' \t '+str(E5HEX[i,0]) +' \t '+ str(E6HEX[i,0]) + ' \t ' + str(E7HEX[i,0])+ ' \t ' +str(E8HEX[i,0]))
    file.write("\n")

    #######################################
    file.write("VARIABLES = \"x\"\"y\"\"z\"\"Rho\"")
    file.write("\n")
    file.write("ZONE STRANDID=3,VARSHARELIST=([1,2,3,4]), SOLUTIONTIME=1, NODES="+str(Nnodes)+", ELEMENTS="+ str(Nelemstet)+", DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON")
 
    for i in range(Nelemstet):
        file.write("\n" + str(E1TET[i,0]) +' \t '+ str(E2TET[i,0]) + ' \t ' + str(E3TET[i,0])+ ' \t ' +str(E4TET[i,0]))
    file.write("\n")
    file.close()

    print("Complete!")
    
