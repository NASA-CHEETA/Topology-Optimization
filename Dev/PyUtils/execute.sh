#!/bin/bash



echo " Initializing Aerodynamic module ..."
echo "               "
mpirun  -n 4 SU2_CFD Fluid/turb_NACA0012_sa.cfg > log.CFD &

echo " Initializing Structural module ..."

time ccx_preCICE -i Solid/wing3d -precice-participant Calculix > log.FEA &

rm -rf flow*
rm -rf res*
rm -rf calculix-*
rm -rf precice-SU2*
rm -rf precice-Calculix-e*
rm -rf surface_flow*
rm -rf precice-SU2_*

wait
