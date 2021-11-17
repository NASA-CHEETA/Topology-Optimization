#!/bin/bash


# echo " Removing old directories ..."

# Remove old files

rm -rf Fluid/surf*
rm -rf Events
rm -rf Probe
rm -rf logs

# echo " Initializing Aerodynamic module ..."
# echo "               "
time mpirun -n 4 SU2_FSI Fluid/inv_flap.cfg >log.CFD 2>&1 &  

# echo " Initializing Structural module ..."
export OMP_NUM_THREADS=2
time ccx_2.15_mkle -i Solid/squareflap -precice-participant Calculix -l 0 >log.FEA 2 >&1 &



wait

# echo " Aeroelastic simulation complete ! "

mkdir Events
cp -r precice-Calc* Events/

mkdir Probe
cp -r precice-Calculix-watchpoint-* Probe/

mkdir logs

cp -r res* Fluid/
cp -r surface* Fluid/

cp -r log.* logs
cp -r history* logs
rm -rf log.*
rm -rf hist*
rm -rf forces_break*
rm -rf flow*
rm -rf res*
rm -rf calculix-*
rm -rf precice-SU2*
rm -rf precice-Calculix-e*
rm -rf surface_flow*
rm -rf precice-SU2_*
rm -rf precice-Cal*
rm -rf spooles.out
rm -rf config_CFD.cfg
rm -rf __pycac*

#echo " Post-processing data ... "

#mpirun -np 14 SU2_SOL2 Fluid/inv_flap.cfg > log.postProcess &

#wait

#echo " All tasks complete ! "

