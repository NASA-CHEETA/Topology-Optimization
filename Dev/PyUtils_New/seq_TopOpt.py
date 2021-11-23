#Script to perform Sequential Aero-elastic Topology Optmization

import os
import shutil  
import subprocess
import glob
import wrtInp
import readwriteData as rwd
import TO_config
import FSI
import meshutil
import ForcSurf

def sTopOpt(Nnodes, Nelems, E, nu, rho, n):

    #Clear all design folders from previous runs
    for file in glob.glob('DSN*'):
        shutil.rmtree(file)
    for file in glob.glob('*temp'):
        shutil.rmtree(file)
    for file in glob.glob('*.inp'):
        os.remove(file)    

    #Get the root directory 
    Root = os.getcwd()

    #Create a temp directory in the first iteration 
    if n == 1:
        os.makedirs("temp")
        os.chdir("temp")
        TMP_DIR = os.getcwd()
        os.chdir(Root)

    #Make TMP_DIR available to everyone 
    TMP_DIR = TMP_DIR 

    #Create a design directory for each design iteration 
    os.makedirs("DSN_00"+str(n))

    #Go to design directory and get its path 
    os.chdir("DSN_00"+str(n))
    DSN_DIR = os.getcwd()

    #With this design directory, create a path for Topology_Opt and FSI_n
    os.makedirs("TopOpt_00"+str(n))
    os.makedirs("FSI_00"+str(n))
    
    #Go into TopOpt Dir and get path
    os.chdir("TopOpt_00"+str(n))
    TOPT_DIR = os.getcwd()

    #Go back to Design Directory and go into FSI_n to get path
    os.chdir(DSN_DIR)
    os.chdir("FSI_00"+str(n))
    FSI_DIR = os.getcwd()

    #Within FSI DIR, create Fluid and Solid sub-directories
    os.makedirs("Fluid")
    os.makedirs("Solid")

    #Get path of Fluid Dir
    os.chdir("Fluid")
    FLUID_DIR = os.getcwd()
    
    os.chdir(FSI_DIR)
    
    #Get path of Solid Dir
    os.chdir("Solid")
    SOLID_DIR = os.getcwd()

    #Go back to root to begin movement of files 
    os.chdir(Root)

    #For the second iteration, take the density_FSI.dat from TMP_DIR and transfer
    # to DSN_DIR
    if n>1:
        os.chdir(TMP_DIR)
        for file in glob.glob("density_FSI.dat"):
            shutil.copy(file, TOPT_DIR)
            shutil.move(file, DSN_DIR)
    #Go back to the root again to begin movement of global files
    os.chdir(Root)

#--------------------------------------------------------------------------------------------------#
#-------------------------------- SUB-DIRECTORY DEFINITIONS COMPLETE ------------------------------#
#--------------------------------------------------------------------------------------------------#

    #Copy all Fluid files from root to FLUID_DIR
    for file in glob.glob("*.su2"):
        shutil.copy(file,FLUID_DIR)

    for file in glob.glob("*.cfg"):
        shutil.copy(file,FLUID_DIR)    

    #Copy all Solid files from root to SOLID_DIR
    for file in glob.glob("*.nam"):
        shutil.copy(file,SOLID_DIR)
    
    #If on first design cycle, get undeformed mesh from Root and create an inp file simultaneously
    if n == 1:
        for file in glob.glob("*.msh"):
            Flag = wrtInp.wrtInput(file,1,E,nu,rho)
            shutil.copy(file, SOLID_DIR)
        for file in glob.glob("*.inp"):
            shutil.move(file,SOLID_DIR)    
    
    #If not on first design cycle, fetch the deformed mesh of previous design cycle from TEMP_DIR
    if n > 1:
        os.chdir(TMP_DIR)
        for file in glob.glob("*.*msh"):
            Flag = wrtInp.wrtInput(file,1,E,nu,rho)
            shutil.move(file,SOLID_DIR)
        for file in glob.glob("*.inp"):
            shutil.move(file,SOLID_DIR)    
        os.chdir(Root)

     #Copy all FSI-Setting files from root to FSI_DIR
    Flag = shutil.copy('precice-config.xml',FSI_DIR)
    Flag = shutil.copy('config.yml',FSI_DIR)
    Flag = shutil.copy('execute.sh',FSI_DIR)


#--------------------------------------------------------------------------------------------------#
#-------------------------------- INITIAL FILE TRANSFER COMPLETE ----------------------------------#
#--------------------------------------------------------------------------------------------------#

    #Move to FSI DIR and perform Aero-elastic Simulation 
    os.chdir(FSI_DIR)

    print("AERO-ELASTIC ITER:", n)
    Flag = FSI.AeroElastic()

    #Read Disp.csv and create deformed mesh 

    for file in glob.glob('Solid/*.msh'):
        Mesh_Name = file  
    Def_Mesh_Name = meshutil.MshUpdt(Mesh_Name, Nnodes, Nelems,n)

    #Move the deformed mesh, Forces and nam files to TopOpt_DIR

    for file in glob.glob('Forces.csv'):
        shutil.move(file,TOPT_DIR)

    os.chdir(SOLID_DIR)
    for file in glob.glob('deform*'):
        shutil.move(file,TOPT_DIR)
    for file in glob.glob('*.nam'):
        shutil.copy(file,TOPT_DIR)
    
    

#--------------------------------------------------------------------------------------------------#
#-------------------------------- AERO-ELASTIC SIMULATION COMPLETE --------------------------------#
#--------------------------------------------------------------------------------------------------#

    #Move into TOPT_DIR for Topology Optimization 
    os.chdir(TOPT_DIR)

    #Create an input file for TopOpt
    for file in glob.glob('*msh'):
        Flag = wrtInp.wrtInput(file,0,E,nu,rho)

    #Create a file for Forces .csv
    Flag = ForcSurf.ForcSurf_to_ccx("Forces.csv", "surfaceForces.nam")


    #Perform Topology Optimization
   # if n == 1:
   #     Flag = TO_config.TopOpt(restartFlag = False)

  #  if n > 1:
        #At this step a density_FSI.dat is read
  #      Flag = TO_config.TopOpt(restartFlag = True)

    #Extract second column from rho.dat for next iter
 #   rhoPhys = rwd.getRhoPhys()
 #   flag   = rwd.writerhoPhys(rhoPhys)

    #Move the density_FSI.dat to tempDir
#    for file in glob.glob("density_FSI.dat"):
#        shutil.move(file, TMP_DIR)

    #Copy the deformed mesh to tempDir for FSI at next design cycle
    for file in glob.glob('*.msh'):
        shutil.copy(file, TMP_DIR)

    #Go back to root at the end 
    os.chdir(Root)


#--------------------------------------------------------------------------------------------------#
#-------------------------------- SEQUENTIAL TOP-OPT DESIGN CYCLE COMPLETE ------------------------#
#--------------------------------------------------------------------------------------------------#





    


