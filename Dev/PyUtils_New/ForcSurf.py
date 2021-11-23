import numpy as np
import sys
import pandas as pd

def ForcSurf_to_ccx(csv, output):
    """
    Takes Fx,Fy,Fz on Nodes in csv file and writes into cload file
    used by ccx inp
    """

    #Extract the node information
    Nodes = pd.read_csv(csv,header = None, sep = ',', skiprows = 1)
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

# Create  global force vectors for all nodes including surface nodes

    FFX =np.zeros((Num_Surf_Nodes,1))
    FFY =np.zeros((Num_Surf_Nodes,1))
    FFZ =np.zeros((Num_Surf_Nodes,1))

    FFX = FFX.reshape(Num_Surf_Nodes,1)
    FFY = FFY.reshape(Num_Surf_Nodes,1)
    FFZ = FFZ.reshape(Num_Surf_Nodes,1)

    file = open(output, "w")
    for i in range(Num_Surf_Nodes):
        file.write( str(int(Nodeid[i,0])) +' ,\t '+ str(1) + ', \t ' + str(FX[i,0]) + ' \n ')
        file.write( str(int(Nodeid[i,0])) +' ,\t '+ str(2) + ', \t ' + str(FY[i,0])+ ' \n ')
        file.write( str(int(Nodeid[i,0])) +' ,\t '+ str(3) + ', \t ' + str(FZ[i,0])+ ' \n ')
    file.write("\n")
    file.close()
    print("Complete!")



if __name__=="__main__":
    
    
    csv= "Forces.csv"   # csv filename given by Prateek
    output= "surfaceNodesopt.nam"   # 
    ForcSurf_to_ccx(csv, output)