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

E   =  2.5E5
nu  =  0.35
rho =  100

# Get rid of old directories if not restarting 
for file in glob.glob('DSN*'):
    shutil.rmtree(file)

for file in glob.glob('temp'):
    shutil.rmtree(file)

for file in glob.glob('*.inp'):
    os.remove(file)


while n < 2:
    n += 1 # get counter for condition to increment if condition is true:

    

    # Create an input deck for CCX for the first iteration only
    if (n ==1):
        for file in glob.glob('*.msh'):
            Flag = wrtInp.wrtInput(file,1,E,nu,rho)
            if ( Flag!=0):
                print(" Error creating an input deck at iteration: ", n)

    # Get the root directory
        Root = os.getcwd()

    # Create a temp directory only at first iteration
    # "temp" holds all the elastic files for next iteration FSI   
        os.makedirs("temp")
        os.chdir("temp")
        Temp_Dir = os.getcwd()
    # Go back to root    
    os.chdir(Root)

    
    # Create a design directory for each design iteration and set its path  
    os.makedirs("DSN_"+str(n))

    # Go into design directoy and get the path
    os.chdir("DSN_"+str(n))
    New_Dir= os.getcwd()

    os.chdir(Root)

    # Copy the FSI-related files to the Design Directory
    
    New_Loc = shutil.copy('precice-config.xml', New_Dir )
    New_Loc = shutil.copy('config.yml', New_Dir )
    New_Loc = shutil.copy('execute.sh', New_Dir )
    
    # Copy the fluid (.su2 + .cfg files) to the DSN directory

    for file in glob.glob('*.cfg'):
        shutil.copy(file, New_Dir)

    for file in glob.glob('*.su2'):
        shutil.copy(file, New_Dir)    

    # Extract the Solid (.inp + .nam + .msh files) to the DSN directory

    for file in glob.glob('*.nam'):
        shutil.copy(file, New_Dir)
    
    # If the first iteration, grab the inp and msh file from root
    if ( n == 1):

        for file in glob.glob('*.inp'):
            shutil.copy(file, New_Dir)
        for file in glob.glob('*.msh'):
            shutil.copy(file, New_Dir) 


    # If not the first iteration, grab the updated inp and msh file from temp

    if (n > 1):
        os.chdir(Temp_Dir)
        for file in glob.glob('*.inp'):
            shutil.move(file, New_Dir)
        for file in glob.glob('deform*'):
            shutil.move(file, New_Dir)
        os.chdir(Root)        
    
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

    print(" FSI Iteration : ", n)

    Flag = FSI.AeroElastic()

    # Read the Disp.csv and create new deformed mesh for topOpt

    # Detect the current mesh in "Solid" for mesh update

    for file in glob.glob('Solid/*.msh'):
        Mesh_Name = file 

    Def_Mesh_Name = meshutil.MshUpdt(Mesh_Name,Nnodes,Nelems,n)


    # Create a sub-directory for FSI Solution

    os.makedirs("FSI_"+str(n))

    FSI_Dir = "FSI_"+str(n)

    # Create a sub-directory for Static Direct Analysis

    os.makedirs("DIRECT_"+str(n))
    DIR_Dir = "DIRECT_"+str(n)


    ## Copy files required for forward analysis

    # Copy the Deformed mesh for static analysis
    for file in glob.glob(Def_Mesh_Name):
        shutil.copy(file, DIR_Dir)

    # Copy the *.nam files for static analysis
    for file in glob.glob('Solid/*.nam'):
        shutil.copy(file, DIR_Dir)     

    # Move the forces.csv for forward analysis
    for file in glob.glob('Forces.csv'):
        shutil.move(file, DIR_Dir)

    # Go to into Forward Analysis Directory

    os.chdir(DIR_Dir)

    # Register path for DIRECT directory 

    DIRECT_Dir = os.getcwd()

    # Create an input deck for CCX static analysis
    for file in glob.glob('*.msh'):
        Flag = wrtInp.wrtInput(file,0,E,nu,rho)
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

    

    for file in glob.glob('Solid'):
        shutil.move(file, FSI_Dir) 

    for file in glob.glob('*.csv'):
        shutil.move(file, FSI_Dir)
    
    # Before moving the deformed mesh file to FSI , copy for next Iteration to temp

   # for file in glob.glob('deform*'):
   #     shutil.copy(file, Temp_Dir )
         
    # Create a sub-directory for Adjoint Sensitivity analysis

    os.makedirs("ADJOINT_"+str(n))

    os.chdir(DIRECT_Dir)
    for file in glob.glob('deform*'):
        shutil.copy(file, Temp_Dir )

     
    # Go to tempDir
    os.chdir(Temp_Dir)

    # Create an input deck for CCX dynamic analysis for next iteration
    for file in glob.glob('*.msh'):
        Flag = wrtInp.wrtInput(file,1,E,nu,rho)
        if ( Flag!=0):
            print(" Error creating an input deck for dynamic analysis at iteration: ", n+1)

    # Go to Root
    os.chdir(Root)







