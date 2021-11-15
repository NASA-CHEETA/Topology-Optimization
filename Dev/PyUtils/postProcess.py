import numpy as np
import pandas as pd
import math
import sys
import displacements
import Forc_read
import Stress_read

TopOpt = False


Nnodes= 10467
Nelems= 45513

#Extract the node information
Nodes = pd.read_csv('Solid/wing_deform.nam',header = None, sep = ',', skiprows = 2, nrows=Nnodes)

#Extract the element information
Element = pd.read_csv('Solid/wing_deform.nam',header = None, sep = ',', skiprows = Nnodes + 7, nrows=Nelems)

if TopOpt:
    # All density-related actions in this block
    #Extract the density information
    Rhos = pd.read_csv('Solid/rhos.dat',header = None, sep = ',', nrows=Nelems)
    RHO = Rhos.to_numpy()
    RHO_FILTERED=RHO[0:Nelems,1]
    RHO_FILTERED=RHO_FILTERED.reshape(Nelems,1)
    RHO_UNFILTERED=RHO[0:Nelems,0]
    RHO_UNFILTERED=RHO_UNFILTERED.reshape(Nelems,1)
    densityinterp = np.zeros(Nnodes)    #To store sum of connecting element rhos at node
    el_to_node_count = np.zeros(Nnodes) #To store number of elements sharing each node
    print("Finding elements belonging to a node...")
    #Look for all elements connecting at a node, and add their rhos
    for elcount, elem in enumerate(Eleid):
        rhoe=RHO_FILTERED[elem-1]

        densityinterp[E1[elem-1]-1] += rhoe
        el_to_node_count[E1[elem-1]-1] += 1

        densityinterp[E2[elem-1]-1] += rhoe
        el_to_node_count[E2[elem-1]-1] += 1

        densityinterp[E3[elem-1]-1] += rhoe
        el_to_node_count[E3[elem-1]-1] += 1

        densityinterp[E4[elem-1]-1] += rhoe
        el_to_node_count[E4[elem-1]-1] += 1


    #Average density at node
    densityinterp_ave = np.divide(densityinterp, el_to_node_count)
    densityinterp_ave = densityinterp_ave.reshape(Nnodes,1)

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
## Extract nodal displacement values
#-----------------------------------------------------------------------------#

U = displacements.displacements(Nnodes)

UX = U[:,1]
UY = U[:,2]
UZ = U[:,3]

UX = UX.reshape(Nnodes,1)
UY = UY.reshape(Nnodes,1)
UZ = UZ.reshape(Nnodes,1)

#-----------------------------------------------------------------------------#
## Extract nodal force values
#-----------------------------------------------------------------------------#

F = Forc_read.Forc_read(Nnodes)
#print(F.shape)
FX = F[:,0]
FY = F[:,1]
FZ = F[:,2]

#print(FY.shape)
#print(FZ.shape)
FX = FX.reshape(Nnodes,1)
FY = FY.reshape(Nnodes,1)
FZ = FZ.reshape(Nnodes,1)
#print(FX.shape)
#print(FY.shape)
#print(FZ.shape)

#-----------------------------------------------------------------------------#
## Extract nodal Sigma values
#-----------------------------------------------------------------------------#

Sigma = Stress_read.Stress_read(Nnodes)




file = open("structure.dat", "w")
file.write("TITLE = \"Visualization of the solid solution\"")
file.write("\n")
if TopOpt:
    file.write("VARIABLES = \"x\"\"y\"\"z\"\"Ux\"\"Uy\"\"Uz\"\"Fx\"\"Fy\"\"Fz\"\"Sigma\"\"Rho\"")
else:
    file.write("VARIABLES = \"x\"\"y\"\"z\"\"Ux\"\"Uy\"\"Uz\"\"Fx\"\"Fy\"\"Fz\"\"Sigma\"")

file.write("\n")
file.write("ZONE STRANDID=2, SOLUTIONTIME=1, NODES= 10467, ELEMENTS= 45513, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON")
for i in range(Nnodes):
    if TopOpt:
        file.write("\n" +"{:.5e}" .format(X[i,0]) +'  \t  '+ "{:.5e}".format(Y[i,0]) + ' \t \t ' + "{:.5e}".format(Z[i,0]) + '\t'+ "{:.5e}".format(UX[i,0])+'\t'+ "{:.5e}".format(UY[i,0])+'\t'+ "{:.5e}".format(UZ[i,0]) + ' \t ' +"{:.5e}".format(FX[i,0]) + ' \t ' +"{:.5e}".format(FY[i,0]) + ' \t ' +"{:.5e}".format(FZ[i,0])+ ' \t ' +"{:.5e}".format(Sigma[i,0])+ ' \t ' +"{:.5e}".format(densityinterp_ave[i,0]))
    else:
        file.write("\n" +"{:.5e}" .format(X[i,0]) +'  \t  '+ "{:.5e}".format(Y[i,0]) + ' \t \t ' + "{:.5e}".format(Z[i,0]) + '\t'+ "{:.5e}".format(UX[i,0])+'\t'+ "{:.5e}".format(UY[i,0])+'\t'+ "{:.5e}".format(UZ[i,0]) + ' \t ' +"{:.5e}".format(FX[i,0]) + ' \t ' +"{:.5e}".format(FY[i,0]) + ' \t ' +"{:.5e}".format(FZ[i,0])+ ' \t ' +"{:.5e}".format(Sigma[i,0]))

for i in range(Nelems):
    file.write("\n" + str(E1[i,0]) +' \t '+ str(E2[i,0]) + ' \t ' + str(E3[i,0])+ ' \t ' +str(E4[i,0]))
file.write("\n")
file.close()
print(" Post processing Complete!")
