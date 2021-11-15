import os
import shutil
import FSI
import meshutil
import subprocess
import glob
import wrtInp

n = 0 # set initial 0 

Nnodes = 7946
Nelems = 31005
while n < 1:
    n += 1 # get counter for condition to increment if condition is true:

    # Create an input deck for CCX for the first iteration only
    if ( n ==1):
        Flag = wrtInp.wrtInput("squareflap.msh")
        if ( Flag!=0):
            print(" Error creating an input deck at iteration: ", n)
        # Get the root directory
        Root = os.getcwd()

        # Make a temporary directory at root
        # "temp" holds all the elastic files for next iteration FSI

        os.makedirs("temp")
        os.chdir("temp")
        Temp_Dir = os.getcwd()
        os.chdir(Root)
    
    # Create a design directory for Optimization and set its path  
    os.makedirs("DSN_"+str(n))
    New_Dir = "DSN_"+str(n)
    
    New_Loc = shutil.copy('precice-config.xml', New_Dir )
    New_Loc = shutil.copy('config.yml', New_Dir )
    New_Loc = shutil.copy('execute.sh', New_Dir )
    
    # Extract the fluid (.su2 files) to the DSN directory

    for file in glob.glob('*.cfg'):
        shutil.copy(file, New_Dir)

    for file in glob.glob('*.su2'):
        shutil.copy(file, New_Dir)    

    # Extract the Solid (.inp, .nam files) to the DSN directory

    for file in glob.glob('*.nam'):
        shutil.copy(file, New_Dir)

    for file in glob.glob('*.inp'):
        shutil.copy(file, New_Dir)

    for file in glob.glob('*.msh'):
        shutil.copy(file, New_Dir)         

    # Enter DSN and re-allocate Fluid and Solid Files   
    os.chdir("DSN_"+str(n))

    # Make Fluid and Solid sub-dirs
    os.makedirs("Fluid")
    os.makedirs("Solid")
    
    # Move the relevant files to Fluid

    for file in glob.glob('*.cfg'):
        shutil.move(file, 'Fluid')

    for file in glob.glob('*.su2'):
        shutil.move(file, 'Fluid') 

    # Move the relevant files to Solid

    for file in glob.glob('*.nam'):
        shutil.move(file, 'Solid')

    for file in glob.glob('*.inp'):
        shutil.move(file, 'Solid')

    for file in glob.glob('*.msh'):
        shutil.move(file, 'Solid')     

    # Perform Forward Analysis in design Directory

    print(" Iteration : ", n)

    Flag = FSI.AeroElastic()

    # Read the Disp.csv and create new deformed mesh for topOpt

    # Detect the current mesh in "Solid"

    for file in glob.glob('Solid/*.msh'):
        Mesh_Name = file 

    Def_Mesh_Name = meshutil.MshUpdt(Mesh_Name,Nnodes,Nelems,n)


    # Create a sub-directory for FSI Solution

    os.makedirs("FSI_"+str(n))

    FSI_Dir = "FSI_"+str(n)

    # Create a sub-directory for Static Direct Analysis

    os.makedirs("DIRECT_"+str(n))
    DIR_Dir = "DIRECT_"+str(n)


    # Get the path again for design directory
    New_Dir = os.getcwd()

    ## Copy files required for forward analysis

    # Copy the Deformed mesh for static analysis
    for file in glob.glob(Def_Mesh_Name):
        shutil.copy(file, DIR_Dir)

    # Copy the *.nam files for static analysis
    for file in glob.glob('Solid/*.nam'):
        shutil.copy(file, DIR_Dir)     

    # CMove the forces.csv for forward analysis
    for file in glob.glob('Forces.csv'):
        shutil.move(file, DIR_Dir)

    # Go to into Forward Analysis Directory

    os.chdir(DIR_Dir)

    # Register path for DIRECT directory 

    DIRECT_Dir = os.getcwd()

    # Create an input deck for CCX static analysis
    for file in glob.glob('*.msh'):
        Flag = wrtInp.wrtInput(file)
        if ( Flag!=0):
            print(" Error creating an input deck for forward analysis at iteration: ", n)

    # Go back to design directory
    os.chdir(New_Dir)        
    
    # All necessary files for foward analysis are acquired

    # Move all FSI files and results inside FSI Dir

    New_Loc = shutil.move('precice-config.xml', FSI_Dir)
    New_Loc = shutil.move('config.yml', FSI_Dir)
    New_Loc = shutil.move('execute.sh', FSI_Dir)
    New_Loc = shutil.move('logs', FSI_Dir)
    New_Loc = shutil.move('Probe', FSI_Dir)
    New_Loc = shutil.move('Events', FSI_Dir)

    for file in glob.glob('Fluid'):
        shutil.move(file, FSI_Dir)

    # Before moving the elastic files to FSI , copy for next Iteration to temp

    for file in glob.glob('Solid/*.nam'):
        shutil.copy(file, Temp_Dir) 

    for file in glob.glob('Solid'):
        shutil.move(file, FSI_Dir) 

    for file in glob.glob('*.csv'):
        shutil.move(file, FSI_Dir)
    
    # Before moving the deformed mesh file to FSI , copy for next Iteration

    for file in glob.glob('*.msh'):
        shutil.copy(file, Temp_Dir )
   

    for file in glob.glob('*.msh'):
        shutil.move(file, FSI_Dir)          

    # Create a sub-directory for Adjoint Sensitivities

    os.makedirs("ADJOINT_"+str(n))

    # Go to tempDir
    os.chdir(Temp_Dir)

    # Create an input deck for CCX static analysis
    for file in glob.glob('*.msh'):
        Flag = wrtInp.wrtInput(file)
        if ( Flag!=0):
            print(" Error creating an input deck for dynamic analysis at iteration: ", n+1)

    # Go to tempDir
    os.chdir(Root)



    

    print("Execution complete!")






