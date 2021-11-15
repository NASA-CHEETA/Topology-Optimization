import numpy as np
import sys
import pandas as pd
import math
import os

#def Forc_read(Nnodes):
#Nnodes = 10467
def Stress_read(Nnodes):
    Str = pd.read_csv('Stress.csv',header = None, sep = ',', skiprows = 3, nrows=Nnodes, engine='python')
    STRESS = Str.to_numpy()
# Determine the number of surface nodes

    SXX =  STRESS[:,1]
    SYY =  STRESS[:,2]
    SZZ =  STRESS[:,3]
    SXY =  STRESS[:,4]
    SYZ =  STRESS[:,5]
    SZX =  STRESS[:,6]


#print(STRESS.shape)
#print(STRESS[:,1])
#print(STRESS[:,2])
#print(STRESS[:,3])

    SXX = SXX.reshape(Nnodes,1)
    SYY = SYY.reshape(Nnodes,1)
    SZZ = SZZ.reshape(Nnodes,1)
    SXY = SXY.reshape(Nnodes,1)
    SYZ = SYZ.reshape(Nnodes,1)
    SZX = SZX.reshape(Nnodes,1)

# Compute von Mises stress for "general" state of stress
    Sigma = np.zeros(Nnodes)
    Sigma = Sigma.reshape(Nnodes,1)
    for i in range(Nnodes):
        Sigma[i,0] = math.sqrt((0.5*((SXX[i,0]-SYY[i,0])**2 + (SYY[i,0]-SZZ[i,0])**2 + (SZZ[i,0]-SXX[i,0])**2 + 3*(SXY[i,0]**2 + SYZ[i,0]**2 + SZX[i,0]**2))))

    os.remove("Stress_read.pyc")
    return Sigma
