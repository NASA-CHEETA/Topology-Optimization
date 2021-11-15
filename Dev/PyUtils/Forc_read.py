import numpy as np
import sys
import pandas as pd

def Forc_read(Nnodes):

    #Extract the node information
    Nodes = pd.read_csv('Forces.csv',header = None, sep = ',', skiprows = 1, nrows=Nnodes)
    LOC = Nodes.to_numpy()
    # Determine the number of surface nodes

    Num_Surf_Nodes = LOC[:,1].size


    Nodeid = LOC[0:Num_Surf_Nodes,0]
    FX = LOC[0:Num_Surf_Nodes,1]
    FY = LOC[0:Num_Surf_Nodes,2]
    FZ = LOC[0:Num_Surf_Nodes,3]


    Nodeid= Nodeid.reshape(Num_Surf_Nodes,1)
    FX = FX.reshape(Num_Surf_Nodes,1)
    FY = FY.reshape(Num_Surf_Nodes,1)
    FZ = FZ.reshape(Num_Surf_Nodes,1)

    if Nodeid[-1]!=Num_Surf_Nodes:
        sys.exit("ERROR : Inconsistent number of surface nodes read from Forces.csv")

# Create  global force vectors for all nodes including surface nodes

    FFX =np.zeros((Nnodes,1))
    FFY =np.zeros((Nnodes,1))
    FFZ =np.zeros((Nnodes,1))

    FFX = FFX.reshape(Nnodes,1)
    FFY = FFY.reshape(Nnodes,1)
    FFZ = FFZ.reshape(Nnodes,1)


    for i in range(Num_Surf_Nodes):
        FFX[i,0] = FX[i,0]
        FFY[i,0] = FY[i,0]
        FFZ[i,0] = FZ[i,0]
    #print(FFX.shape)
    #print(FFY.shape)
    #print(FFZ.shape)


    F = np.concatenate((FFX,FFY,FFZ),axis = 1)
    return F
