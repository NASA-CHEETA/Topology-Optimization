#!/usr/bin/env python

## \file serial_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 6.0.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function, division, absolute_import
import sys
from TestCase import TestCase

def main():
    '''This program runs SU2 and ensures that the output matches specified values. 
       This will be used to do checks when code is pushed to github 
       to make sure nothing is broken. '''

    test_list = []
    
    #########################
    ## Compressible Euler ###
    #########################

    # Channel
    channel           = TestCase('channel')
    channel.cfg_dir   = "euler/channel"
    channel.cfg_file  = "inv_channel_RK.cfg"
    channel.test_iter = 10
    channel.test_vals = [-2.454049, 3.065639, -0.200679, 0.036298] #last 4 columns
    channel.su2_exec  = "SU2_CFD"
    channel.timeout   = 1600
    channel.tol       = 0.00001
    test_list.append(channel)

    # NACA0012 
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.047448, -3.538057, 0.338691, 0.023131] #last 4 columns
    naca0012.su2_exec  = "SU2_CFD"
    naca0012.timeout   = 1600
    naca0012.tol       = 0.00001
    test_list.append(naca0012)

    # Supersonic wedge 
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-0.804690, 4.936631, -0.251188, 0.044255] #last 4 columns
    wedge.su2_exec  = "SU2_CFD"
    wedge.timeout   = 1600
    wedge.tol       = 0.00001
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-13.393130, -12.928941, 0.282557, 0.012706] #last 4 columns
    oneram6.su2_exec  = "SU2_CFD"
    oneram6.timeout   = 9600
    oneram6.tol       = 0.00001
    test_list.append(oneram6)
    
    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 100
    fixedCL_naca0012.test_vals = [-2.582897, 2.820443, 0.295158, 0.019324] #last 4 columns
    fixedCL_naca0012.su2_exec  = "SU2_CFD"
    fixedCL_naca0012.timeout   = 1600
    fixedCL_naca0012.tol       = 0.00001
    test_list.append(fixedCL_naca0012)
    
    # Polar sweep of the inviscid NACA0012
    polar_naca0012           = TestCase('polar_naca0012')
    polar_naca0012.cfg_dir   = "polar/naca0012"
    polar_naca0012.cfg_file  = "inv_NACA0012.cfg"
    polar_naca0012.polar     = True
    polar_naca0012.test_iter = 10
    polar_naca0012.test_vals = [-1.319488, 4.112397, 0.011954, 0.009584] #last 4 columns
    polar_naca0012.su2_exec  = "compute_polar.py -n 1 -i 11"
    polar_naca0012.timeout   = 1600
    polar_naca0012.tol       = 0.00001
    test_list.append(polar_naca0012)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 20
    flatplate.test_vals = [-4.680911, 0.781118, -0.135956, 0.022978] #last 4 columns
    flatplate.su2_exec  = "SU2_CFD"
    flatplate.timeout   = 1600
    flatplate.tol       = 0.00001
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-6.765426, -1.297422, 0.019496, 0.310082] #last 4 columns
    cylinder.su2_exec  = "SU2_CFD"
    cylinder.timeout   = 1600
    cylinder.tol       = 0.00001
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-6.850123, -1.388088, -0.056090, 108.140177] #last 4 columns
    cylinder_lowmach.su2_exec  = "SU2_CFD"
    cylinder_lowmach.timeout   = 1600
    cylinder_lowmach.tol       = 0.00001
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-12.272126, -3.335311, 0.000001, 2.351005] #last 4 columns
    poiseuille.su2_exec  = "SU2_CFD"
    poiseuille.timeout   = 1600
    poiseuille.tol       = 0.00001
    test_list.append(poiseuille)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-2.000469, -5.228296, 0.820188, 0.052004] #last 4 columns
    rae2822_sa.su2_exec  = "SU2_CFD"
    rae2822_sa.timeout   = 1600
    rae2822_sa.tol       = 0.00001
    test_list.append(rae2822_sa)
    
    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-0.510826, 4.909749, 0.825022, 0.052675] #last 4 columns
    rae2822_sst.su2_exec  = "SU2_CFD"
    rae2822_sst.timeout   = 1600
    rae2822_sst.tol       = 0.00001
    test_list.append(rae2822_sst)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.158341, -6.737135, -0.176245, 0.057446] #last 4 columns
    turb_flatplate.su2_exec  = "SU2_CFD"
    turb_flatplate.timeout   = 1600
    turb_flatplate.tol       = 0.00001
    test_list.append(turb_flatplate)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.327509, -6.563372, 0.230438, 0.155815]#last 4 columns
    turb_oneram6.su2_exec  = "SU2_CFD"
    turb_oneram6.timeout   = 3200
    turb_oneram6.tol       = 0.00001
    test_list.append(turb_oneram6)

    # NACA0012 (SA, FUN3D results for finest grid: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 10
    turb_naca0012_sa.test_vals = [-12.000762, -9.145363, 1.070528, 0.019417] #last 4 columns
    turb_naca0012_sa.su2_exec  = "SU2_CFD"
    turb_naca0012_sa.timeout   = 3200
    turb_naca0012_sa.tol       = 0.00001
    test_list.append(turb_naca0012_sa)
   
    # NACA0012 (SA, FUN3D results for finest grid: CL=1.0983, CD=0.01242) with binary restart
    turb_naca0012_sa_bin           = TestCase('turb_naca0012_sa_bin')
    turb_naca0012_sa_bin.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa_bin.cfg_file  = "turb_NACA0012_sa_binary.cfg"
    turb_naca0012_sa_bin.test_iter = 10
    turb_naca0012_sa_bin.test_vals = [-12.000899, -9.145363, 1.070528, 0.019417] #last 4 columns
    turb_naca0012_sa_bin.su2_exec  = "SU2_CFD"
    turb_naca0012_sa_bin.timeout   = 3200
    turb_naca0012_sa_bin.tol       = 0.00001
    test_list.append(turb_naca0012_sa_bin)
 
    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-15.039675, -7.219913, 1.059622, 0.019138] #last 4 columns
    turb_naca0012_sst.su2_exec  = "SU2_CFD"
    turb_naca0012_sst.timeout   = 3200
    turb_naca0012_sst.tol       = 0.00001
    test_list.append(turb_naca0012_sst)

    # PROPELLER 
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.378804, -8.130405, 0.000051, 0.055764] #last 4 columns
    propeller.su2_exec  = "SU2_CFD"
    propeller.timeout   = 3200
    propeller.tol       = 0.00001
    test_list.append(propeller)

    #################################
    ## Compressible RANS Restart  ###
    #################################

    # NACA0012 SST Multigrid restart
    turb_naca0012_sst_restart_mg           = TestCase('turb_naca0012_sst_restart_mg')
    turb_naca0012_sst_restart_mg.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_restart_mg.cfg_file  = "turb_NACA0012_sst_multigrid_restart.cfg"
    turb_naca0012_sst_restart_mg.test_iter = 50
    turb_naca0012_sst_restart_mg.ntest_vals = 5
    turb_naca0012_sst_restart_mg.test_vals = [-6.459444, -4.595710, 1.201844, -0.007146, 0.080517] #last 5 columns
    turb_naca0012_sst_restart_mg.su2_exec  = "SU2_CFD"
    turb_naca0012_sst_restart_mg.timeout   = 3200
    turb_naca0012_sst_restart_mg.tol       = 0.000001
    test_list.append(turb_naca0012_sst_restart_mg)


    #############################
    ### Incompressible Euler  ###
    #############################

    # NACA0012 Hydrofoil
    inc_euler_naca0012           = TestCase('inc_euler_naca0012')
    inc_euler_naca0012.cfg_dir   = "incomp_euler/naca0012"
    inc_euler_naca0012.cfg_file  = "incomp_NACA0012.cfg"
    inc_euler_naca0012.test_iter = 20
    inc_euler_naca0012.test_vals = [-3.590333,-3.151515,0.960192,0.013260] #last 4 columns
    inc_euler_naca0012.su2_exec  = "SU2_CFD"
    inc_euler_naca0012.timeout   = 1600
    inc_euler_naca0012.tol       = 0.00001
    test_list.append(inc_euler_naca0012)

    #############################
    ### Incompressible N-S    ###
    #############################

    # Laminar cylinder
    inc_lam_cylinder          = TestCase('inc_lam_cylinder')
    inc_lam_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_lam_cylinder.cfg_file  = "incomp_cylinder.cfg"
    inc_lam_cylinder.test_iter = 10
    inc_lam_cylinder.test_vals = [-3.778049,-3.380095,-0.004749,6.215563] #last 4 columns
    inc_lam_cylinder.su2_exec  = "SU2_CFD"
    inc_lam_cylinder.timeout   = 1600
    inc_lam_cylinder.tol       = 0.00001
    test_list.append(inc_lam_cylinder)

    ############################
    ### Incompressible RANS  ###
    ############################

    # NACA0012
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.708169, -11.010140, -0.000002, 0.207687] #last 4 columns
    inc_turb_naca0012.su2_exec  = "SU2_CFD"
    inc_turb_naca0012.timeout   = 1600
    inc_turb_naca0012.tol       = 0.00001
    test_list.append(inc_turb_naca0012)

    #########################
    ###    Transition     ###
    #########################

    # Schubauer-Klebanoff Natural Transition
    schubauer_klebanoff_transition              = TestCase('Schubauer_Klebanoff')
    schubauer_klebanoff_transition.cfg_dir      = "transition/Schubauer_Klebanoff"
    schubauer_klebanoff_transition.cfg_file     = "transitional_BC_model_ConfigFile.cfg"
    schubauer_klebanoff_transition.test_iter    = 10
    schubauer_klebanoff_transition.test_vals    = [-8.295514, -14.505574, 0.000073, 0.005669] #last 4 columns
    schubauer_klebanoff_transition.su2_exec     = "SU2_CFD"
    schubauer_klebanoff_transition.timeout      = 1600
    schubauer_klebanoff_transition.tol          = 0.00001
    test_list.append(schubauer_klebanoff_transition)

    #####################################
    ### Cont. adj. compressible Euler ###
    #####################################

    # Inviscid NACA0012
    contadj_naca0012           = TestCase('contadj_naca0012')
    contadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    contadj_naca0012.cfg_file  = "inv_NACA0012.cfg"
    contadj_naca0012.test_iter = 5
    contadj_naca0012.test_vals = [-9.787554, -15.192510, 0.300920, 0.019552] #last 4 columns
    contadj_naca0012.su2_exec  = "SU2_CFD"
    contadj_naca0012.timeout   = 1600
    contadj_naca0012.tol       = 0.001
    test_list.append(contadj_naca0012)

    # Inviscid ONERA M6
    contadj_oneram6           = TestCase('contadj_oneram6')
    contadj_oneram6.cfg_dir   = "cont_adj_euler/oneram6"
    contadj_oneram6.cfg_file  = "inv_ONERAM6.cfg"
    contadj_oneram6.test_iter = 10
    contadj_oneram6.test_vals = [-12.133352, -12.707213, 6.8590e-01, 0.007594] #last 4 columns
    contadj_oneram6.su2_exec  = "SU2_CFD"
    contadj_oneram6.timeout   = 1600
    contadj_oneram6.tol       = 0.00001
    test_list.append(contadj_oneram6)

    # Inviscid WEDGE: tests averaged outflow total pressure adjoint
    contadj_wedge             = TestCase('contadj_wedge')
    contadj_wedge.cfg_dir   = "cont_adj_euler/wedge"
    contadj_wedge.cfg_file  = "inv_wedge_ROE.cfg"
    contadj_wedge.test_iter = 10
    contadj_wedge.test_vals = [2.856008, -2.767216, 1.0029e+06, 1.3024e-13] #last 4 columns
    contadj_wedge.su2_exec  = "SU2_CFD"
    contadj_wedge.timeout   = 1600
    contadj_wedge.tol       = 0.00001
    test_list.append(contadj_wedge)
    
    # Inviscid fixed CL NACA0012
    contadj_fixedCL_naca0012           = TestCase('contadj_fixedcl_naca0012')
    contadj_fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    contadj_fixedCL_naca0012.cfg_file  = "inv_NACA0012_ContAdj.cfg"
    contadj_fixedCL_naca0012.test_iter = 100
    contadj_fixedCL_naca0012.test_vals = [0.341038, -5.166613, 0.265510, -0.000322] #last 4 columns
    contadj_fixedCL_naca0012.su2_exec  = "SU2_CFD"
    contadj_fixedCL_naca0012.timeout   = 1600
    contadj_fixedCL_naca0012.tol       = 0.00001
    test_list.append(contadj_fixedCL_naca0012)

    ###################################
    ### Cont. adj. compressible N-S ###
    ###################################

    # Adjoint laminar cylinder
    contadj_ns_cylinder           = TestCase('contadj_ns_cylinder')
    contadj_ns_cylinder.cfg_dir   = "cont_adj_navierstokes/cylinder"
    contadj_ns_cylinder.cfg_file  = "lam_cylinder.cfg"
    contadj_ns_cylinder.test_iter = 20
    contadj_ns_cylinder.test_vals = [ -3.661071, -9.126563, 2.056700, -0.000000] #last 4 columns
    contadj_ns_cylinder.su2_exec  = "SU2_CFD"
    contadj_ns_cylinder.timeout   = 1600
    contadj_ns_cylinder.tol       = 0.00001
    test_list.append(contadj_ns_cylinder)

    # Adjoint laminar naca0012 subsonic
    contadj_ns_naca0012_sub           = TestCase('contadj_ns_naca0012_sub')
    contadj_ns_naca0012_sub.cfg_dir   = "cont_adj_navierstokes/naca0012_sub"
    contadj_ns_naca0012_sub.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_sub.test_iter = 20
    contadj_ns_naca0012_sub.test_vals = [-2.743268, -8.215193, 0.518810, 0.001210] #last 4 columns
    contadj_ns_naca0012_sub.su2_exec  = "SU2_CFD"
    contadj_ns_naca0012_sub.timeout   = 1600
    contadj_ns_naca0012_sub.tol       = 0.00001
    test_list.append(contadj_ns_naca0012_sub)
    
    # Adjoint laminar naca0012 transonic
    contadj_ns_naca0012_trans           = TestCase('contadj_ns_naca0012_trans')
    contadj_ns_naca0012_trans.cfg_dir   = "cont_adj_navierstokes/naca0012_trans"
    contadj_ns_naca0012_trans.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_trans.test_iter = 20
    contadj_ns_naca0012_trans.test_vals = [-1.039664, -6.575019, 1.772300, 0.012495] #last 4 columns
    contadj_ns_naca0012_trans.su2_exec  = "SU2_CFD"
    contadj_ns_naca0012_trans.timeout   = 1600
    contadj_ns_naca0012_trans.tol       = 0.00001
    test_list.append(contadj_ns_naca0012_trans)

    #######################################################
    ### Cont. adj. compressible RANS (frozen viscosity) ###
    #######################################################

    # Adjoint turbulent NACA0012
    contadj_rans_naca0012           = TestCase('contadj_rans_naca0012')
    contadj_rans_naca0012.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012.cfg_file  = "turb_nasa.cfg"
    contadj_rans_naca0012.test_iter = 20
    contadj_rans_naca0012.test_vals = [ -0.794162, -5.761722, 19.214000, -0.000000] #last 4 columns
    contadj_rans_naca0012.su2_exec  = "SU2_CFD"
    contadj_rans_naca0012.timeout   = 1600
    contadj_rans_naca0012.tol       = 0.00001
    test_list.append(contadj_rans_naca0012)
   
    # Adjoint turbulent NACA0012 with binary restarts
    contadj_rans_naca0012_bin           = TestCase('contadj_rans_naca0012_bin')
    contadj_rans_naca0012_bin.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012_bin.cfg_file  = "turb_nasa_binary.cfg"
    contadj_rans_naca0012_bin.test_iter = 20
    contadj_rans_naca0012_bin.test_vals = [-0.794169, -5.761671, 19.214000, -0.000000] #last 4 columns
    contadj_rans_naca0012_bin.su2_exec  = "SU2_CFD"
    contadj_rans_naca0012_bin.timeout   = 1600
    contadj_rans_naca0012_bin.tol       = 0.00001
    test_list.append(contadj_rans_naca0012_bin)
 
    # Adjoint turbulent RAE2822
    contadj_rans_rae2822           = TestCase('contadj_rans_rae2822')
    contadj_rans_rae2822.cfg_dir   = "cont_adj_rans/rae2822"
    contadj_rans_rae2822.cfg_file  = "turb_SA_RAE2822.cfg"
    contadj_rans_rae2822.test_iter = 20
    contadj_rans_rae2822.test_vals = [-5.369688, -10.872209, -0.212470, 0.005448] #last 4 columns
    contadj_rans_rae2822.su2_exec  = "SU2_CFD"
    contadj_rans_rae2822.timeout   = 1600
    contadj_rans_rae2822.tol       = 0.00001
    test_list.append(contadj_rans_rae2822)

    #######################################
    ### Cont. adj. incompressible Euler ###
    #######################################

    # Adjoint Incompressible Inviscid NACA0012
    contadj_incomp_NACA0012           = TestCase('contadj_incomp_NACA0012')
    contadj_incomp_NACA0012.cfg_dir   = "cont_adj_incomp_euler/naca0012"
    contadj_incomp_NACA0012.cfg_file  = "incomp_NACA0012.cfg"
    contadj_incomp_NACA0012.test_iter = 5
    contadj_incomp_NACA0012.test_vals = [-11.980272, -12.146779, 1.9399, 0.000000] #last 4 columns
    contadj_incomp_NACA0012.su2_exec  = "SU2_CFD"
    contadj_incomp_NACA0012.timeout   = 1600
    contadj_incomp_NACA0012.tol       = 0.00001
    test_list.append(contadj_incomp_NACA0012)

    #####################################
    ### Cont. adj. incompressible N-S ###
    #####################################

    # Adjoint Incompressible Viscous Cylinder
    contadj_incomp_cylinder           = TestCase('contadj_incomp_cylinder')
    contadj_incomp_cylinder.cfg_dir   = "cont_adj_incomp_navierstokes/cylinder"
    contadj_incomp_cylinder.cfg_file  = "lam_incomp_cylinder.cfg"
    contadj_incomp_cylinder.test_iter = 25
    contadj_incomp_cylinder.test_vals = [-5.720871, -6.995398, 2.932100, 0.000000] #last 4 columns
    contadj_incomp_cylinder.su2_exec  = "SU2_CFD"
    contadj_incomp_cylinder.timeout   = 1600
    contadj_incomp_cylinder.tol       = 0.00001
    test_list.append(contadj_incomp_cylinder)

    ######################################
    ### Harmonic Balance               ###
    ######################################
    
    # Description of the regression test
    harmonic_balance           = TestCase('harmonic_balance')
    harmonic_balance.cfg_dir   = "harmonic_balance"
    harmonic_balance.cfg_file  = "HB.cfg"
    harmonic_balance.test_iter = 25
    harmonic_balance.test_vals = [-1.569573, 3.941896, 0.008780, 0.079775] #last 4 columns
    harmonic_balance.su2_exec  = "SU2_CFD"
    harmonic_balance.timeout   = 1600
    harmonic_balance.tol       = 0.00001
    test_list.append(harmonic_balance)

    # Turbulent pitching NACA 64a010 airfoil
    hb_rans_preconditioning           = TestCase('hb_rans_preconditioning')
    hb_rans_preconditioning.cfg_dir   = "harmonic_balance/hb_rans_preconditioning"
    hb_rans_preconditioning.cfg_file  = "davis.cfg"
    hb_rans_preconditioning.test_iter = 25
    hb_rans_preconditioning.test_vals = [-1.900976,-5.876722, 0.007761, 0.125934] #last 4 columns
    hb_rans_preconditioning.su2_exec  = "SU2_CFD"
    hb_rans_preconditioning.timeout   = 1600
    hb_rans_preconditioning.tol       = 0.00001
    test_list.append(hb_rans_preconditioning)

    ######################################
    ### Moving Wall                    ###
    ######################################
    
    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-5.627934, -0.164470, 0.051972, 2.547034] #last 4 columns
    cavity.su2_exec  = "SU2_CFD"
    cavity.timeout   = 1600
    cavity.tol       = 0.00001
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals = [-7.709662, -2.274900, 1.418422, 1.734206] #last 4 columns
    spinning_cylinder.su2_exec  = "SU2_CFD"
    spinning_cylinder.timeout   = 1600
    spinning_cylinder.tol       = 0.00001
    test_list.append(spinning_cylinder)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-1.166421, 0.076804, 1.398548, 2.197047] #last 4 columns
    square_cylinder.su2_exec  = "SU2_CFD"
    square_cylinder.timeout   = 1600
    square_cylinder.tol       = 0.00001
    square_cylinder.unsteady  = True
    test_list.append(square_cylinder)

    # Gust
    sine_gust           = TestCase('sine_gust')
    sine_gust.cfg_dir   = "gust"
    sine_gust.cfg_file  = "inv_gust_NACA0012.cfg"
    sine_gust.test_iter = 5
    sine_gust.test_vals = [-1.977531, 3.481790, -0.006222, -0.001342] #last 4 columns
    sine_gust.su2_exec  = "SU2_CFD"
    sine_gust.timeout   = 1600
    sine_gust.tol       = 0.00001
    sine_gust.unsteady  = True
    test_list.append(sine_gust)

    # Aeroelastic
    aeroelastic         = TestCase('aeroelastic')
    aeroelastic.cfg_dir   = "aeroelastic"
    aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    aeroelastic.test_iter = 2
    aeroelastic.test_vals = [0.077106, 0.036449, -1.684916e-03, -1.131735e-04] #last 4 columns
    aeroelastic.su2_exec  = "SU2_CFD"
    aeroelastic.timeout   = 1600
    aeroelastic.tol       = 0.000001
    aeroelastic.unsteady  = True
    test_list.append(aeroelastic) 

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714721, -5.883008, -0.214968, 0.023783] #last 4 columns
    ddes_flatplate.su2_exec  = "SU2_CFD"
    ddes_flatplate.timeout   = 1600
    ddes_flatplate.tol       = 0.00001
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)    

    ######################################
    ### NICFD                          ###
    ######################################

    # ls89_sa
    ls89_sa           = TestCase('ls89_sa')
    ls89_sa.cfg_dir   = "nicf/LS89"
    ls89_sa.cfg_file  = "turb_SA_PR.cfg"
    ls89_sa.test_iter = 20
    ls89_sa.test_vals = [-5.050658, -13.393681, 0.174911, 0.430700] #last 4 columns
    ls89_sa.su2_exec  = "SU2_CFD"
    ls89_sa.timeout   = 1600
    ls89_sa.tol       = 0.00001
    test_list.append(ls89_sa)

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 20
    edge_VW.test_vals = [-0.720910, 5.481122, -0.000853, 0.000000] #last 4 columns
    edge_VW.su2_exec  = "SU2_CFD"
    edge_VW.timeout   = 1600
    edge_VW.tol       = 0.00001
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR               
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 20
    edge_PPR.test_vals = [-1.688033, 4.505202, 0.000916, 0.000000] #last 4 columns
    edge_PPR.su2_exec  = "SU2_CFD"
    edge_PPR.timeout   = 1600
    edge_PPR.tol       = 0.00001
    test_list.append(edge_PPR)
    
    
    ######################################
    ### turbomachinery                 ###
    ######################################
    
    # Jones APU Turbocharger
    Jones_tc           = TestCase('jones_turbocharger')
    Jones_tc.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc.cfg_file  = "Jones.cfg"
    Jones_tc.test_iter = 5
    Jones_tc.test_vals = [-5.362540, 0.406893, 82.840290, 0.882243] #last 4 columns
    Jones_tc.su2_exec  = "SU2_CFD"
    Jones_tc.timeout   = 1600
    Jones_tc.tol       = 0.00001
    test_list.append(Jones_tc)

	# Jones APU Turbocharger restart
    Jones_tc_rst           = TestCase('jones_turbocharger_restart')
    Jones_tc_rst.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc_rst.cfg_file  = "Jones_rst.cfg"
    Jones_tc_rst.test_iter = 5
    Jones_tc_rst.test_vals = [-4.265376, -1.422420, 82.253590, 2.776443] #last 4 columns
    Jones_tc_rst.su2_exec  = "SU2_CFD"
    Jones_tc_rst.timeout   = 1600
    Jones_tc_rst.tol       = 0.00001
    test_list.append(Jones_tc_rst)

    # 2D axial stage
    axial_stage2D           = TestCase('axial_stage2D')
    axial_stage2D.cfg_dir   = "turbomachinery/axial_stage_2D"
    axial_stage2D.cfg_file  = "Axial_stage2D.cfg"
    axial_stage2D.test_iter = 20
    axial_stage2D.test_vals = [-1.837052, 5.787278, 73.679940, 0.888933] #last 4 columns
    axial_stage2D.su2_exec  = "SU2_CFD"
    axial_stage2D.timeout   = 1600
    axial_stage2D.tol       = 0.00001
    test_list.append(axial_stage2D)
    
    # 2D transonic stator
    transonic_stator           = TestCase('transonic_stator')
    transonic_stator.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator.cfg_file  = "transonic_stator.cfg"
    transonic_stator.test_iter = 20
    transonic_stator.test_vals = [-1.111211, 6.120416, 67.101490, 0.070037] #last 4 columns
    transonic_stator.su2_exec  = "SU2_CFD"
    transonic_stator.timeout   = 1600
    transonic_stator.tol       = 0.00001
    test_list.append(transonic_stator)
    
    # 2D transonic stator restart
    transonic_stator_rst           = TestCase('transonic_stator_restart')
    transonic_stator_rst.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator_rst.cfg_file  = "transonic_stator_rst.cfg"
    transonic_stator_rst.test_iter = 20
    transonic_stator_rst.test_vals = [-0.487105, 4.460638, 6.456525, 0.004009] #last 4 columns
    transonic_stator_rst.su2_exec  = "SU2_CFD"
    transonic_stator_rst.timeout   = 1600
    transonic_stator_rst.tol       = 0.00001
    test_list.append(transonic_stator_rst)


    ######################################
    ### Sliding Mesh                   ###
    ######################################

    # Uniform flow
    uniform_flow         = TestCase('uniform_flow')
    uniform_flow.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow.cfg_file  = "uniform_NN.cfg"
    uniform_flow.test_iter = 50
    uniform_flow.test_vals = [-0.368836, 5.156090, 0.000000, 0.000000] #last 4 columns
    uniform_flow.su2_exec  = "SU2_CFD"
    uniform_flow.timeout   = 1600
    uniform_flow.tol       = 0.000001
    uniform_flow.unsteady  = True
    test_list.append(uniform_flow) 

    # Channel_2D
    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 4
    channel_2D.test_vals = [-1.656896, 4.263158, 0.000000, 0.000000] #last 4 columns
    channel_2D.su2_exec  = "SU2_CFD"
    channel_2D.timeout   = 100
    channel_2D.tol       = 0.00001
    channel_2D.unsteady  = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 1
    channel_3D.test_vals = [-1.899946, 4.032882, 0.000000, 0.000000] #last 4 columns
    channel_3D.su2_exec  = "SU2_CFD"
    channel_3D.timeout   = 1600
    channel_3D.tol       = 0.00001
    channel_3D.unsteady  = True
    test_list.append(channel_3D)

    # Pipe
    pipe           = TestCase('pipe')
    pipe.cfg_dir   = "sliding_interface/pipe"
    pipe.cfg_file  = "pipe_NN.cfg"
    pipe.test_iter = 2
    pipe.test_vals = [-3.503708, 3.194241, 0.000000, 0.000000] #last 4 columns
    pipe.su2_exec  = "SU2_CFD"
    pipe.timeout   = 1600
    pipe.tol       = 0.00001
    pipe.unsteady  = True
    test_list.append(pipe)

    # Rotating cylinders
    rotating_cylinders           = TestCase('rotating_cylinders')
    rotating_cylinders.cfg_dir   = "sliding_interface/rotating_cylinders"
    rotating_cylinders.cfg_file  = "rot_cylinders_WA.cfg"
    rotating_cylinders.test_iter = 3
    rotating_cylinders.test_vals = [-1.253451, 4.531335, 0.000000, 0.000000] #last 4 columns
    rotating_cylinders.su2_exec  = "SU2_CFD"
    rotating_cylinders.timeout   = 1600
    rotating_cylinders.tol       = 0.00001
    rotating_cylinders.unsteady  = True
    test_list.append(rotating_cylinders)

    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [-1.122916, 4.606248, 0.000000, 0.000000] #last 4 columns
    supersonic_vortex_shedding.su2_exec  = "SU2_CFD"
    supersonic_vortex_shedding.timeout   = 1600
    supersonic_vortex_shedding.tol       = 0.00001
    supersonic_vortex_shedding.unsteady  = True
    test_list.append(supersonic_vortex_shedding)

    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [-2.138791, 1.639495, -0.000831, 0.117497] #last 4 columns
    bars_SST_2D.su2_exec  = "SU2_CFD"
    bars_SST_2D.timeout   = 1600
    bars_SST_2D.tol       = 0.00001
    test_list.append(bars_SST_2D)

    ##########################
    ### FEA - FSI          ###
    ##########################

    # Static beam, 3d
    statbeam3d           = TestCase('statbeam3d')
    statbeam3d.cfg_dir   = "fea_fsi/StatBeam_3d"
    statbeam3d.cfg_file  = "configBeam_3d.cfg"
    statbeam3d.test_iter = 0
    statbeam3d.test_vals = [-8.498274, -8.230638, -8.123824, 6.4095e+04] #last 4 columns
    statbeam3d.su2_exec  = "SU2_CFD"
    statbeam3d.timeout   = 1600
    statbeam3d.tol       = 0.00001
    test_list.append(statbeam3d)

    # Dynamic beam, 2d
    dynbeam2d           = TestCase('dynbeam2d')
    dynbeam2d.cfg_dir   = "fea_fsi/DynBeam_2d"
    dynbeam2d.cfg_file  = "configBeam_2d.cfg"
    dynbeam2d.test_iter = 6
    dynbeam2d.test_vals = [-9.420640, -5.365872, -12.430382, 6.5210e+04] #last 4 columns
    dynbeam2d.su2_exec  = "SU2_CFD"
    dynbeam2d.timeout   = 1600
    dynbeam2d.tol       = 0.00001
    test_list.append(dynbeam2d)

    # FSI, 2d
    fsi2d           = TestCase('fsi2d')
    fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    fsi2d.cfg_file  = "configFSI_2D.cfg"
    fsi2d.test_iter = 4
    fsi2d.test_vals = [2.000000, 0.500000, -7.780236, -1.142100] #last 4 columns
    fsi2d.su2_exec  = "SU2_CFD"
    fsi2d.timeout   = 1600
    fsi2d.tol       = 0.00001
    test_list.append(fsi2d)    
   

    ######################################
    ### RUN TESTS                      ###
    ######################################  

    pass_list = [ test.run_test() for test in test_list ]

    
    ######################################
    ### RUN SU2_GEO TESTS              ###
    ######################################
    
    # NACA0012
    naca0012_geo           = TestCase('naca0012_geo')
    naca0012_geo.cfg_dir   = "optimization_euler/steady_naca0012"
    naca0012_geo.cfg_file  = "inv_NACA0012_adv.cfg"
    naca0012_geo.test_vals = [0.99975, 62.0455, 0.120041, 0.00208202] #chord, LE radius, ToC, Alpha
    naca0012_geo.su2_exec  = "SU2_GEO"
    naca0012_geo.timeout   = 1600
    naca0012_geo.tol       = 0.00001
    pass_list.append(naca0012_geo.run_geo())
    test_list.append(naca0012_geo)

    ######################################
    ### RUN SU2_DEF TESTS              ###
    ######################################
    
    # Inviscid NACA0012 (triangles)
    naca0012_def            = TestCase('naca0012_def')
    naca0012_def.cfg_dir   = "deformation/naca0012"
    naca0012_def.cfg_file  = "def_NACA0012.cfg"
    naca0012_def.test_iter = 10
    naca0012_def.test_vals = [4.604230e-03] #residual
    naca0012_def.su2_exec  = "SU2_DEF"
    naca0012_def.timeout   = 1600
    naca0012_def.tol       = 1e-08
    
    pass_list.append(naca0012_def.run_def())
    test_list.append(naca0012_def)
    
    # RAE2822 (mixed tris + quads)
    rae2822_def            = TestCase('rae2822_def')
    rae2822_def.cfg_dir   = "deformation/rae2822"
    rae2822_def.cfg_file  = "def_RAE2822.cfg"
    rae2822_def.test_iter = 10
    rae2822_def.test_vals = [8.348610e-09] #residual
    rae2822_def.su2_exec  = "SU2_DEF"
    rae2822_def.timeout   = 1600
    rae2822_def.tol       = 1e-13
    
    pass_list.append(rae2822_def.run_def())
    test_list.append(rae2822_def)
    
    # Turb NACA4412 (quads, wall distance)
    naca4412_def            = TestCase('naca4412_def')
    naca4412_def.cfg_dir   = "deformation/naca4412"
    naca4412_def.cfg_file  = "def_NACA4412.cfg"
    naca4412_def.test_iter = 10
    naca4412_def.test_vals = [8.855370e-13] #residual
    naca4412_def.su2_exec  = "SU2_DEF"
    naca4412_def.timeout   = 1600
    naca4412_def.tol       = 1e-12
    
    pass_list.append(naca4412_def.run_def())
    test_list.append(naca4412_def)

    # Brick of tets (inverse volume)
    brick_tets_def            = TestCase('brick_tets_def')
    brick_tets_def.cfg_dir   = "deformation/brick_tets"
    brick_tets_def.cfg_file  = "def_brick_tets.cfg"
    brick_tets_def.test_iter = 10
    brick_tets_def.test_vals = [8.973010e-04] #residual
    brick_tets_def.su2_exec  = "SU2_DEF"
    brick_tets_def.timeout   = 1600
    brick_tets_def.tol       = 1e-09
    
    pass_list.append(brick_tets_def.run_def())
    test_list.append(brick_tets_def)

    # Brick of isotropic hexas (inverse volume)
    brick_hex_def           = TestCase('brick_hex_def')
    brick_hex_def.cfg_dir   = "deformation/brick_hex"
    brick_hex_def.cfg_file  = "def_brick_hex.cfg"
    brick_hex_def.test_iter = 10
    brick_hex_def.test_vals = [2.082100e-04] #residual
    brick_hex_def.su2_exec  = "SU2_DEF"
    brick_hex_def.timeout   = 1600
    brick_hex_def.tol       = 1e-09
    
    pass_list.append(brick_hex_def.run_def())
    test_list.append(brick_hex_def)

    # Brick with a pyramid layer (inverse volume)
    brick_pyra_def           = TestCase('brick_pyra_def')
    brick_pyra_def.cfg_dir   = "deformation/brick_pyra"
    brick_pyra_def.cfg_file  = "def_brick_pyra.cfg"
    brick_pyra_def.test_iter = 10
    brick_pyra_def.test_vals = [1.994770e-03] #residual
    brick_pyra_def.su2_exec  = "SU2_DEF"
    brick_pyra_def.timeout   = 1600
    brick_pyra_def.tol       = 1e-08
    
    pass_list.append(brick_pyra_def.run_def())
    test_list.append(brick_pyra_def)

    # Brick of isotropic prisms (inverse volume)
    brick_prism_def           = TestCase('brick_prism_def')
    brick_prism_def.cfg_dir   = "deformation/brick_prism"
    brick_prism_def.cfg_file  = "def_brick_prism.cfg"
    brick_prism_def.test_iter = 10
    brick_prism_def.test_vals = [5.684440e-03] #residual
    brick_prism_def.su2_exec  = "SU2_DEF"
    brick_prism_def.timeout   = 1600
    brick_prism_def.tol       = 1e-08
    
    pass_list.append(brick_prism_def.run_def())
    test_list.append(brick_prism_def)
    
    # Brick of prisms with high aspect ratio cells near the wall (wall distance)
    brick_prism_rans_def           = TestCase('brick_prism_rans_def')
    brick_prism_rans_def.cfg_dir   = "deformation/brick_prism_rans"
    brick_prism_rans_def.cfg_file  = "def_brick_prism_rans.cfg"
    brick_prism_rans_def.test_iter = 10
    brick_prism_rans_def.test_vals = [2.762780e-07] #residual
    brick_prism_rans_def.su2_exec  = "SU2_DEF"
    brick_prism_rans_def.timeout   = 1600
    brick_prism_rans_def.tol       = 1e-12
    
    pass_list.append(brick_prism_rans_def.run_def())
    test_list.append(brick_prism_rans_def)
    
    # Brick of hexas with high aspect ratio cells near the wall (inverse volume)
    brick_hex_rans_def           = TestCase('brick_hex_rans_def')
    brick_hex_rans_def.cfg_dir   = "deformation/brick_hex_rans"
    brick_hex_rans_def.cfg_file  = "def_brick_hex_rans.cfg"
    brick_hex_rans_def.test_iter = 10
    brick_hex_rans_def.test_vals = [2.260750e-07] #residual
    brick_hex_rans_def.su2_exec  = "SU2_DEF"
    brick_hex_rans_def.timeout   = 1600
    brick_hex_rans_def.tol       = 1e-12
    
    pass_list.append(brick_hex_rans_def.run_def())
    test_list.append(brick_hex_rans_def)

    # Cylindrical FFD test
    cylinder_ffd_def           = TestCase('cylinder_ffd_def')
    cylinder_ffd_def.cfg_dir   = "deformation/cylindrical_ffd"
    cylinder_ffd_def.cfg_file  = "def_cylindrical.cfg"
    cylinder_ffd_def.test_iter = 10
    cylinder_ffd_def.test_vals = [4.030020e-04] #residual
    cylinder_ffd_def.su2_exec  = "SU2_DEF"
    cylinder_ffd_def.timeout   = 1600
    cylinder_ffd_def.tol       = 1e-09

    pass_list.append(cylinder_ffd_def.run_def())
    test_list.append(cylinder_ffd_def)

    # Spherical FFD test
    sphere_ffd_def           = TestCase('sphere_ffd_def')
    sphere_ffd_def.cfg_dir   = "deformation/spherical_ffd"
    sphere_ffd_def.cfg_file  = "def_spherical.cfg"
    sphere_ffd_def.test_iter = 10
    sphere_ffd_def.test_vals = [3.733370e-03] #residual
    sphere_ffd_def.su2_exec  = "SU2_DEF"
    sphere_ffd_def.timeout   = 1600
    sphere_ffd_def.tol       = 1e-08

    pass_list.append(sphere_ffd_def.run_def())
    test_list.append(sphere_ffd_def)
    
    # Spherical FFD test using BSplines
    sphere_ffd_def_bspline           = TestCase('sphere_ffd_def_bspline')
    sphere_ffd_def_bspline.cfg_dir   = "deformation/spherical_ffd"
    sphere_ffd_def_bspline.cfg_file  = "def_spherical_bspline.cfg"
    sphere_ffd_def_bspline.test_iter = 10
    sphere_ffd_def_bspline.test_vals = [2.112070e-03] #residual
    sphere_ffd_def_bspline.su2_exec  = "SU2_DEF"
    sphere_ffd_def_bspline.timeout   = 1600
    sphere_ffd_def_bspline.tol       = 1e-08

    pass_list.append(sphere_ffd_def_bspline.run_def())
    test_list.append(sphere_ffd_def_bspline)
   
    ######################################
    ### RUN PYTHON TESTS               ###
    ###################################### 
    
    # test continuous_adjoint.py
    contadj_euler_py = TestCase('contadj_euler_py')
    contadj_euler_py.cfg_dir = "cont_adj_euler/naca0012"
    contadj_euler_py.cfg_file  = "inv_NACA0012.cfg"
    contadj_euler_py.test_iter = 10
    contadj_euler_py.su2_exec  = "continuous_adjoint.py"
    contadj_euler_py.timeout   = 1600
    contadj_euler_py.reference_file = "of_grad_cd.dat.ref"
    contadj_euler_py.test_file = "of_grad_cd.dat"
    pass_list.append(contadj_euler_py.run_filediff())
    test_list.append(contadj_euler_py)

#    # test finite_difference.py
#    findiff_euler_py = TestCase('findiff_euler_py')
#    findiff_euler_py.cfg_dir = "cont_adj_euler/naca0012"
#    findiff_euler_py.cfg_file  = "inv_NACA0012_FD.cfg"
#    findiff_euler_py.test_iter = 10
#    findiff_euler_py.su2_exec  = "finite_differences.py"
#    findiff_euler_py.timeout   = 1600
#    findiff_euler_py.reference_file = "of_grad_findiff.dat.ref"
#    findiff_euler_py.test_file = "FINDIFF/of_grad_findiff.dat"
#    pass_list.append(findiff_euler_py.run_filediff())
#    test_list.append(findiff_euler_py)
    
    # test shape_optimization.py
    shape_opt_euler_py           = TestCase('shape_opt_euler_py')
    shape_opt_euler_py.cfg_dir   = "optimization_euler/steady_naca0012"
    shape_opt_euler_py.cfg_file  = "inv_NACA0012_adv.cfg"
    shape_opt_euler_py.test_iter = 1
    shape_opt_euler_py.test_vals = [1, 1, 2.134974E-05, 3.829535E-03] #last 4 columns
    shape_opt_euler_py.su2_exec  = "shape_optimization.py -g CONTINUOUS_ADJOINT -f"
    shape_opt_euler_py.timeout   = 1600
    shape_opt_euler_py.tol       = 0.00001
    pass_list.append(shape_opt_euler_py.run_opt())
    test_list.append(shape_opt_euler_py)

    # Multiple functionals with the continuous adjoint
    contadj_multi_py            = TestCase('contadj_multi_py')
    contadj_multi_py.cfg_dir    = "cont_adj_euler/wedge"
    contadj_multi_py.cfg_file   = "inv_wedge_ROE_multiobj.cfg"
    contadj_multi_py.test_iter  = 10
    contadj_multi_py.su2_exec   = "continuous_adjoint.py"
    contadj_multi_py.timeout    = 1600
    contadj_multi_py.reference_file = "of_grad_combo.dat.ref"
    contadj_multi_py.test_file  = "of_grad_combo.dat"
    pass_list.append(contadj_multi_py.run_filediff())
    test_list.append(contadj_multi_py)
    
    # Optimization with multiple objectives, with gradients evaluated individually
    # the difference in gradient value relative to combined case 
    # is due to lack of solution file for the adjoint and small number of iterations
    opt_multiobj_py            = TestCase('opt_multiobj_py')
    opt_multiobj_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
    opt_multiobj_py.cfg_file   = "inv_wedge_ROE_multiobj.cfg"
    opt_multiobj_py.test_iter  = 1
    opt_multiobj_py.test_vals = [1, 1, 1.084701E+02, 3.799222E+00] #last 4 columns
    opt_multiobj_py.su2_exec   = "shape_optimization.py -g CONTINUOUS_ADJOINT -f"
    opt_multiobj_py.timeout    = 1600
    opt_multiobj_py.tol       = 0.00001
    pass_list.append(opt_multiobj_py.run_opt())
    test_list.append(opt_multiobj_py)

    # test optimization, with multiple objectives and gradient evaluated as 'combo' 
    opt_multiobjcombo_py            = TestCase('opt_multiobjcombo_py')
    opt_multiobjcombo_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
    opt_multiobjcombo_py.cfg_file   = "inv_wedge_ROE_multiobj_combo.cfg"
    opt_multiobjcombo_py.test_iter  = 1
    opt_multiobjcombo_py.test_vals = [1, 1, 1.084701E+02, 3.789322E+00] #last 4 columns
    opt_multiobjcombo_py.su2_exec   = "shape_optimization.py -g CONTINUOUS_ADJOINT -f"
    opt_multiobjcombo_py.timeout    = 1600
    opt_multiobjcombo_py.tol       = 0.00001
    pass_list.append(opt_multiobjcombo_py.run_opt())
    test_list.append(opt_multiobjcombo_py)

    # test optimization, with multiple objectives evaluated on a single surface
    opt_multiobj1surf_py            = TestCase('opt_multiobj1surf_py')
    opt_multiobj1surf_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
    opt_multiobj1surf_py.cfg_file   = "inv_wedge_ROE_multiobj_1surf.cfg"
    opt_multiobj1surf_py.test_iter  = 1
    opt_multiobj1surf_py.test_vals = [1, 1, 3.083034E+01, 3.789380E+00] #last 4 columns
    opt_multiobj1surf_py.su2_exec   = "shape_optimization.py -g CONTINUOUS_ADJOINT  -f"
    opt_multiobj1surf_py.timeout    = 1600
    opt_multiobj1surf_py.tol       = 0.00001
    pass_list.append(opt_multiobj1surf_py.run_opt())
    test_list.append(opt_multiobj1surf_py)

    # test optimization, with a single objective evaluated on multiple surfaces
    opt_2surf1obj_py            = TestCase('opt_2surf1obj_py')
    opt_2surf1obj_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
    opt_2surf1obj_py.cfg_file   = "inv_wedge_ROE_2surf_1obj.cfg"
    opt_2surf1obj_py.test_iter  = 1    
    opt_2surf1obj_py.test_vals = [1, 1, 2.005657E+00, 1.779489E-07] #last 4 columns
    opt_2surf1obj_py.su2_exec   = "shape_optimization.py -g CONTINUOUS_ADJOINT  -f"
    opt_2surf1obj_py.timeout    = 1600
    opt_2surf1obj_py.tol       = 0.00001
    pass_list.append(opt_2surf1obj_py.run_opt())
    test_list.append(opt_2surf1obj_py)

    ##########################
    ###   Python wrapper   ###
    ##########################
    
    # NACA0012 
    pywrapper_naca0012           = TestCase('pywrapper_naca0012')
    pywrapper_naca0012.cfg_dir   = "euler/naca0012"
    pywrapper_naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    pywrapper_naca0012.test_iter = 20
    pywrapper_naca0012.test_vals = [-4.047448, -3.538057, 0.338691, 0.023131] #last 4 columns
    pywrapper_naca0012.su2_exec  = "SU2_CFD.py -f"
    pywrapper_naca0012.timeout   = 1600
    pywrapper_naca0012.tol       = 0.00001
    test_list.append(pywrapper_naca0012)
    pass_list.append(pywrapper_naca0012.run_test())

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    pywrapper_turb_naca0012_sst           = TestCase('pywrapper_turb_naca0012_sst')
    pywrapper_turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    pywrapper_turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    pywrapper_turb_naca0012_sst.test_iter = 10
    pywrapper_turb_naca0012_sst.test_vals = [-15.039675, -7.219913, 1.059622, 0.019138] #last 4 columns
    pywrapper_turb_naca0012_sst.su2_exec  = "SU2_CFD.py -f"
    pywrapper_turb_naca0012_sst.timeout   = 3200
    pywrapper_turb_naca0012_sst.tol       = 0.00001
    test_list.append(pywrapper_turb_naca0012_sst)
    pass_list.append(pywrapper_turb_naca0012_sst.run_test())

    # Square cylinder
    pywrapper_square_cylinder           = TestCase('pywrapper_square_cylinder')
    pywrapper_square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    pywrapper_square_cylinder.cfg_file  = "turb_square.cfg"
    pywrapper_square_cylinder.test_iter = 3
    pywrapper_square_cylinder.test_vals = [-1.166421, 0.076804, 1.398548, 2.197047] #last 4 columns
    pywrapper_square_cylinder.su2_exec  = "SU2_CFD.py -f"
    pywrapper_square_cylinder.timeout   = 1600
    pywrapper_square_cylinder.tol       = 0.00001
    pywrapper_square_cylinder.unsteady  = True
    test_list.append(pywrapper_square_cylinder)
    pass_list.append(pywrapper_square_cylinder.run_test())

    # Aeroelastic
    pywrapper_aeroelastic         = TestCase('pywrapper_aeroelastic')
    pywrapper_aeroelastic.cfg_dir   = "aeroelastic"
    pywrapper_aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    pywrapper_aeroelastic.test_iter = 2
    pywrapper_aeroelastic.test_vals = [0.077106, 0.036449, -1.684916e-03, -1.131735e-04] #last 4 columns
    pywrapper_aeroelastic.su2_exec  = "SU2_CFD.py -f"
    pywrapper_aeroelastic.timeout   = 1600
    pywrapper_aeroelastic.tol       = 0.000001
    pywrapper_aeroelastic.unsteady  = True
    test_list.append(pywrapper_aeroelastic)
    pass_list.append(pywrapper_aeroelastic.run_test())

    # FSI, 2d
    pywrapper_fsi2d           = TestCase('pywrapper_fsi2d')
    pywrapper_fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    pywrapper_fsi2d.cfg_file  = "configFSI_2D.cfg"
    pywrapper_fsi2d.test_iter = 4
    pywrapper_fsi2d.test_vals = [2.000000, 0.500000, -7.780236, -1.142100] #last 4 columns
    pywrapper_fsi2d.su2_exec  = "SU2_CFD.py --nZone 2 --fsi True -f"
    pywrapper_fsi2d.timeout   = 1600
    pywrapper_fsi2d.tol       = 0.00001
    test_list.append(pywrapper_fsi2d)
    pass_list.append(pywrapper_fsi2d.run_test())

    # Unsteady CHT
    pywrapper_unsteadyCHT               = TestCase('pywrapper_unsteadyCHT')
    pywrapper_unsteadyCHT.cfg_dir       = "py_wrapper/flatPlate_unsteady_CHT"
    pywrapper_unsteadyCHT.cfg_file      = "unsteady_CHT_FlatPlate_Conf.cfg"
    pywrapper_unsteadyCHT.test_iter     = 5
    pywrapper_unsteadyCHT.test_vals     = [-1.598116, 2.263309, 0.001077, 0.145818] #last 4 columns
    pywrapper_unsteadyCHT.su2_exec      = "python launch_unsteady_CHT_FlatPlate.py -f"
    pywrapper_unsteadyCHT.timeout       = 1600
    pywrapper_unsteadyCHT.tol           = 0.00001
    pywrapper_unsteadyCHT.unsteady      = True
    test_list.append(pywrapper_unsteadyCHT)
    pass_list.append(pywrapper_unsteadyCHT.run_test())

    # Rigid motion
    pywrapper_rigidMotion               = TestCase('pywrapper_rigidMotion')
    pywrapper_rigidMotion.cfg_dir       = "py_wrapper/flatPlate_rigidMotion"
    pywrapper_rigidMotion.cfg_file      = "flatPlate_rigidMotion_Conf.cfg"
    pywrapper_rigidMotion.test_iter     = 5
    pywrapper_rigidMotion.test_vals     = [-1.598116, 2.259671, -0.040621, 0.144134] #last 4 columns
    pywrapper_rigidMotion.su2_exec      = "python launch_flatPlate_rigidMotion.py -f"
    pywrapper_rigidMotion.timeout       = 1600
    pywrapper_rigidMotion.tol           = 0.00001
    pywrapper_rigidMotion.unsteady      = True
    test_list.append(pywrapper_rigidMotion)
    pass_list.append(pywrapper_rigidMotion.run_test())
    
    # Tests summary
    print('==================================================================')
    print('Summary of the serial tests')
    print('python version:', sys.version)
    for i, test in enumerate(test_list):
        if (pass_list[i]):
            print('  passed - %s'%test.tag)
        else:
            print('* FAILED - %s'%test.tag)
    
    if all(pass_list):
        sys.exit(0)
    else:
        sys.exit(1)
    # done

if __name__ == '__main__':
    main()
