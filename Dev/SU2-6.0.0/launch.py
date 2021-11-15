import os
import subprocess

for filename in os.listdir(os.getcwd()):
    print filename
    proc = subprocess.Popen(["time ccx_preCICE -i Solid/flap -precice-participant Calculix ", filename])
    proc.wait()
