/*!
 * \file solution_direct_mean.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon
 * \version 6.0.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/solver_structure.hpp"

CEulerSolver::CEulerSolver(void) : CSolver() {
  
  /*--- Basic array initialization ---*/
  
  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  CoPx_Inv = NULL; CoPy_Inv = NULL; CoPz_Inv = NULL;

  CD_Mnt = NULL; CL_Mnt = NULL; CSF_Mnt = NULL;  CEff_Mnt = NULL;
  CMx_Mnt = NULL; CMy_Mnt = NULL; CMz_Mnt = NULL;
  CFx_Mnt = NULL; CFy_Mnt = NULL; CFz_Mnt = NULL;
  CoPx_Mnt = NULL; CoPy_Mnt = NULL; CoPz_Mnt = NULL;

  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  ForceInviscid = NULL; MomentInviscid = NULL;
  ForceMomentum = NULL; MomentMomentum = NULL;

  /*--- Surface based array initialization ---*/
  
  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL_Mnt = NULL; Surface_CD_Mnt = NULL; Surface_CSF_Mnt = NULL; Surface_CEff_Mnt = NULL;
  Surface_CFx_Mnt = NULL; Surface_CFy_Mnt = NULL; Surface_CFz_Mnt = NULL;
  Surface_CMx_Mnt = NULL; Surface_CMy_Mnt = NULL; Surface_CMz_Mnt = NULL;

  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  /*--- Rotorcraft simulation array initialization ---*/
  
  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;
  
  CMerit_Mnt = NULL;  CT_Mnt = NULL;  CQ_Mnt = NULL;

  /*--- Supersonic simulation array initialization ---*/
  
  CEquivArea_Inv = NULL;
  CNearFieldOF_Inv = NULL;
  
  /*--- Engine simulation array initialization ---*/
  
  Inflow_MassFlow = NULL;   Inflow_Pressure = NULL;
  Inflow_Mach = NULL;       Inflow_Area = NULL;
  Exhaust_Pressure = NULL;  Exhaust_Temperature = NULL;
  Exhaust_MassFlow = NULL;  Exhaust_Area = NULL;
  
  /*--- Numerical methods array initialization ---*/
  
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;

  DonorPrimVar = NULL; DonorGlobalIndex = NULL;
  ActDisk_DeltaP = NULL; ActDisk_DeltaT = NULL;

  Smatrix = NULL; Cvector = NULL;
 
  Secondary = NULL; Secondary_i = NULL; Secondary_j = NULL;

  /*--- Fixed CL mode initialization (cauchy criteria) ---*/
  
  Cauchy_Value   = 0;
  Cauchy_Func    = 0;
  Old_Func       = 0;
  New_Func       = 0;
  Cauchy_Counter = 0;
  Cauchy_Serie = NULL;
  
  AoA_FD_Change = false;

  FluidModel   = NULL;
  
  SlidingState     = NULL;
  SlidingStateNodes = NULL;

  /*--- Initialize quantities for the average process for internal flow ---*/
  
  AverageVelocity            = NULL;
  AverageTurboVelocity       = NULL;
  ExtAverageTurboVelocity    = NULL;
  AverageFlux                = NULL;
  SpanTotalFlux              = NULL;
  AveragePressure            = NULL;
  RadialEquilibriumPressure  = NULL;
  ExtAveragePressure         = NULL;
  AverageDensity             = NULL;
  ExtAverageDensity          = NULL;
  AverageNu                  = NULL;
  AverageKine                = NULL;
  AverageOmega               = NULL;
  ExtAverageNu               = NULL;
  ExtAverageKine             = NULL;
  ExtAverageOmega            = NULL;


  /*--- Initialize primitive quantities for turboperformace ---*/

  DensityIn                     = NULL;
  PressureIn                    = NULL;
  TurboVelocityIn               = NULL;
  DensityOut                    = NULL;
  PressureOut                   = NULL;
  TurboVelocityOut              = NULL;
  KineIn                        = NULL;
  OmegaIn                       = NULL;
  NuIn                          = NULL;
  KineOut                       = NULL;
  OmegaOut                      = NULL;
  NuOut                         = NULL;


  /*--- Initialize quantities for Giles BC---*/

  CkInflow                      = NULL;
  CkOutflow1                    = NULL;
  CkOutflow2                    = NULL;
 
}

CEulerSolver::CEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  
  unsigned long iPoint, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  su2double StaticEnergy, Density, Velocity2, Pressure, Temperature;
  unsigned short nZone = geometry->GetnZone();
  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  bool roe_turkel = (config->GetKind_Upwind_Flow() == TURKEL);
  bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS));
  unsigned short direct_diff = config->GetDirectDiff();
  int Unst_RestartIter;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  bool low_mach_prec = config->Low_Mach_Preconditioning();

  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  string filename_ = config->GetSolution_FlowFileName();

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone);

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/

    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Read and store the restart metadata. ---*/

    Read_SU2_Restart_Metadata(geometry, config, false, filename_);

  }

  /*--- Array initialization ---*/

  /*--- Basic array initialization ---*/

  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  CoPx_Inv = NULL; CoPy_Inv = NULL; CoPz_Inv = NULL;

  CD_Mnt= NULL; CL_Mnt= NULL; CSF_Mnt= NULL; CEff_Mnt= NULL;
  CMx_Mnt= NULL;   CMy_Mnt= NULL;   CMz_Mnt= NULL;
  CFx_Mnt= NULL;   CFy_Mnt= NULL;   CFz_Mnt= NULL;
  CoPx_Mnt= NULL;   CoPy_Mnt= NULL;   CoPz_Mnt= NULL;

  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  ForceInviscid = NULL; MomentInviscid = NULL;
  ForceMomentum = NULL;  MomentMomentum = NULL;

  /*--- Surface based array initialization ---*/
  
  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL_Mnt= NULL; Surface_CD_Mnt= NULL; Surface_CSF_Mnt= NULL; Surface_CEff_Mnt= NULL;
  Surface_CFx_Mnt= NULL;   Surface_CFy_Mnt= NULL;   Surface_CFz_Mnt= NULL;
  Surface_CMx_Mnt= NULL;   Surface_CMy_Mnt= NULL;   Surface_CMz_Mnt = NULL;

  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  /*--- Rotorcraft simulation array initialization ---*/

  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;

  CMerit_Mnt = NULL; CT_Mnt = NULL; CQ_Mnt = NULL;

  /*--- Supersonic simulation array initialization ---*/
  
  CEquivArea_Inv = NULL;
  CNearFieldOF_Inv = NULL;
  
  /*--- Engine simulation array initialization ---*/
  
  Inflow_MassFlow = NULL;   Inflow_Pressure = NULL;
  Inflow_Mach = NULL;       Inflow_Area = NULL;
  Exhaust_Pressure = NULL;  Exhaust_Temperature = NULL;
  Exhaust_MassFlow = NULL;  Exhaust_Area = NULL;
  
  /*--- Numerical methods array initialization ---*/
  
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;
  DonorPrimVar = NULL; DonorGlobalIndex = NULL;
  ActDisk_DeltaP = NULL; ActDisk_DeltaT = NULL;

  Smatrix = NULL; Cvector = NULL;

  Secondary=NULL; Secondary_i=NULL; Secondary_j=NULL;

  /*--- Fixed CL mode initialization (cauchy criteria) ---*/

  Cauchy_Value = 0;
  Cauchy_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  Cauchy_Counter = 0;
  Cauchy_Serie = NULL;
  
  AoA_FD_Change = false;

  FluidModel = NULL;
  
  /*--- Initialize quantities for the average process for internal flow ---*/

  AverageVelocity                   = NULL;
  AverageTurboVelocity              = NULL;
  ExtAverageTurboVelocity           = NULL;
  AverageFlux                       = NULL;
  SpanTotalFlux                     = NULL;
  AveragePressure                   = NULL;
  RadialEquilibriumPressure         = NULL;
  ExtAveragePressure                = NULL;
  AverageDensity                    = NULL;
  ExtAverageDensity                 = NULL;
  AverageNu                         = NULL;
  AverageKine                       = NULL;
  AverageOmega                      = NULL;
  ExtAverageNu                      = NULL;
  ExtAverageKine                    = NULL;
  ExtAverageOmega                   = NULL;


  /*--- Initialize primitive quantities for turboperformace ---*/

  DensityIn                     = NULL;
  PressureIn                    = NULL;
  TurboVelocityIn               = NULL;
  DensityOut                    = NULL;
  PressureOut                   = NULL;
  TurboVelocityOut              = NULL;
  KineIn                        = NULL;
  OmegaIn                       = NULL;
  NuIn                          = NULL;
  KineOut                       = NULL;
  OmegaOut                      = NULL;
  NuOut                         = NULL;


  /*--- Initialize quantities for Giles BC---*/

  CkInflow                      = NULL;
  CkOutflow1                    = NULL;
  CkOutflow2                    = NULL;

  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure
   Compressible flow, primitive variables (T, vx, vy, vz, P, rho, h, c, lamMu, EddyMu, ThCond, Cp).
   ---*/
  
  nDim = geometry->GetnDim();

  nVar = nDim+2;
  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
  nSecondaryVar = 2; nSecondaryVarGrad = 2;

  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nPrimVarGrad;
  
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
 
  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex = new unsigned long[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) 
    nVertex[iMarker] = geometry->nVertex[iMarker];
 
  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  
  SetNondimensionalization(geometry, config, iMesh);
  
  /*--- Allocate the node variables ---*/
  
  node = new CVariable*[nPoint];
  
  /*--- Define some auxiliary vectors related to the residual ---*/
  
  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Residual_i    = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
  
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Define some auxiliary vectors related to the primitive solution ---*/
  
  Primitive   = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the Secondary solution ---*/
  
  Secondary   = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar]   = 0.0;
  Secondary_i = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_i[iVar] = 0.0;
  Secondary_j = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }
  
  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/
  
  if (roe_turkel || low_mach_prec) {
    LowMach_Precontioner = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new su2double[nVar];
  }
  
  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  }
  
  else {
    if (rank == MASTER_NODE) cout << "Explicit scheme. No Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Define some auxiliary vectors for computing flow variable
   gradients by least squares, S matrix := inv(R)*traspose(inv(R)),
   c vector := transpose(WA)*(Wb) ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    Cvector = new su2double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      Cvector[iVar] = new su2double [nDim];
    
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  CharacPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }
  
  /*--- Store the value of the primitive variables + 2 turb variables at the boundaries,
   used for IO with a donor cell ---*/
  
  DonorPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      if (rans) {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar+2];
        for (iVar = 0; iVar < nPrimVar + 2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
      else {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
        for (iVar = 0; iVar < nPrimVar ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
    }
  }
  
  /*--- Store the value of the characteristic primitive variables index at the boundaries ---*/
  
  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorGlobalIndex[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Delta P at the Actuator Disk ---*/
  
  ActDisk_DeltaP = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ActDisk_DeltaP[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      ActDisk_DeltaP[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Delta T at the Actuator Disk ---*/
  
  ActDisk_DeltaT = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ActDisk_DeltaT[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      ActDisk_DeltaT[iMarker][iVertex] = 0;
    }
  }

  /*--- Store the value of the Total Pressure at the inlet BC ---*/

  Inlet_Ttotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
      Inlet_Ttotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        Inlet_Ttotal[iMarker][iVertex] = 0;
      }
    } else {
      Inlet_Ttotal[iMarker] = NULL;
    }
  }

  /*--- Store the value of the Total Temperature at the inlet BC ---*/

  Inlet_Ptotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
      Inlet_Ptotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        Inlet_Ptotal[iMarker][iVertex] = 0;
      }
    } else {
      Inlet_Ptotal[iMarker] = NULL;
    }
  }

  /*--- Store the value of the Flow direction at the inlet BC ---*/

  Inlet_FlowDir = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
      Inlet_FlowDir[iMarker] = new su2double* [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        Inlet_FlowDir[iMarker][iVertex] = new su2double [nDim];
        for (iDim = 0; iDim < nDim; iDim++) {
          Inlet_FlowDir[iMarker][iVertex][iDim] = 0;
        }
      }
    } else {
      Inlet_FlowDir[iMarker] = NULL;
    }
  }

  /*--- Set up inlet profiles, if necessary ---*/

  SetInlet(config);

  /*--- Force definition and coefficient arrays for all of the markers ---*/
  
  CPressure = new su2double* [nMarker];
  CPressureTarget = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new su2double [geometry->nVertex[iMarker]];
    CPressureTarget[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressure[iMarker][iVertex] = 0.0;
      CPressureTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Non-dimensional coefficients ---*/
  
  ForceInviscid     = new su2double[nDim];
  MomentInviscid    = new su2double[3];
  CD_Inv         = new su2double[nMarker];
  CL_Inv         = new su2double[nMarker];
  CSF_Inv        = new su2double[nMarker];
  CMx_Inv           = new su2double[nMarker];
  CMy_Inv           = new su2double[nMarker];
  CMz_Inv           = new su2double[nMarker];
  CEff_Inv          = new su2double[nMarker];
  CFx_Inv           = new su2double[nMarker];
  CFy_Inv           = new su2double[nMarker];
  CFz_Inv           = new su2double[nMarker];
  CoPx_Inv           = new su2double[nMarker];
  CoPy_Inv           = new su2double[nMarker];
  CoPz_Inv           = new su2double[nMarker];
 
  ForceMomentum     = new su2double[nDim];
  MomentMomentum    = new su2double[3];
  CD_Mnt        = new su2double[nMarker];
  CL_Mnt        = new su2double[nMarker];
  CSF_Mnt       = new su2double[nMarker];
  CMx_Mnt          = new su2double[nMarker];
  CMy_Mnt          = new su2double[nMarker];
  CMz_Mnt          = new su2double[nMarker];
  CEff_Mnt         = new su2double[nMarker];
  CFx_Mnt          = new su2double[nMarker];
  CFy_Mnt          = new su2double[nMarker];
  CFz_Mnt          = new su2double[nMarker];
  CoPx_Mnt          = new su2double[nMarker];
  CoPy_Mnt          = new su2double[nMarker];
  CoPz_Mnt          = new su2double[nMarker];

  Surface_CL_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Inv     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL_Mnt     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Mnt     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Mnt= new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Mnt        = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF         = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff           = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz            = new su2double[config->GetnMarker_Monitoring()];

  /*--- Rotorcraft coefficients ---*/
  
  CT_Inv           = new su2double[nMarker];
  CQ_Inv           = new su2double[nMarker];
  CMerit_Inv       = new su2double[nMarker];
  
  CT_Mnt           = new su2double[nMarker];
  CQ_Mnt           = new su2double[nMarker];
  CMerit_Mnt       = new su2double[nMarker];

  /*--- Supersonic coefficients ---*/
  
  CEquivArea_Inv   = new su2double[nMarker];
  CNearFieldOF_Inv = new su2double[nMarker];
  
  /*--- Engine simulation ---*/
  
  Inflow_MassFlow     = new su2double[nMarker];
  Inflow_Pressure     = new su2double[nMarker];
  Inflow_Mach         = new su2double[nMarker];
  Inflow_Area         = new su2double[nMarker];
  
  Exhaust_MassFlow    = new su2double[nMarker];
  Exhaust_Pressure    = new su2double[nMarker];
  Exhaust_Temperature = new su2double[nMarker];
  Exhaust_Area        = new su2double[nMarker];
  
  /*--- Init total coefficients ---*/
  
  Total_CD      = 0.0;    Total_CL           = 0.0;    Total_CSF          = 0.0;
  Total_CMx     = 0.0;    Total_CMy          = 0.0;    Total_CMz          = 0.0;
  Total_CoPx    = 0.0;    Total_CoPy         = 0.0;    Total_CoPz         = 0.0;
  Total_CEff    = 0.0;    Total_CEquivArea   = 0.0;    Total_CNearFieldOF = 0.0;
  Total_CFx     = 0.0;    Total_CFy          = 0.0;    Total_CFz          = 0.0;
  Total_CT      = 0.0;    Total_CQ           = 0.0;    Total_CMerit       = 0.0;
  Total_MaxHeat = 0.0;    Total_Heat         = 0.0;    Total_ComboObj     = 0.0;
  Total_CpDiff  = 0.0;    Total_HeatFluxDiff = 0.0;    Total_Custom_ObjFunc=0.0;
  Total_NetCThrust = 0.0; Total_NetCThrust_Prev = 0.0; Total_BCThrust_Prev = 0.0;
  Total_Power = 0.0;      AoA_Prev           = 0.0;
  Total_CL_Prev = 0.0;    Total_CD_Prev      = 0.0;
  Total_CMx_Prev      = 0.0; Total_CMy_Prev      = 0.0; Total_CMz_Prev      = 0.0;
  Total_AeroCD = 0.0;  Total_SolidCD = 0.0;   Total_IDR   = 0.0;    Total_IDC   = 0.0;

  /*--- Read farfield conditions ---*/
  
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Mach_Inf        = config->GetMach();
  
  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  
  switch(direct_diff) {
    case NO_DERIVATIVE:
      /*--- Default ---*/
      break;
    case D_DENSITY:
      SU2_TYPE::SetDerivative(Density_Inf, 1.0);
      break;
    case D_PRESSURE:
      SU2_TYPE::SetDerivative(Pressure_Inf, 1.0);
      break;
    case D_TEMPERATURE:
      SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
      break;
    case D_MACH: case D_AOA:
    case D_SIDESLIP: case D_REYNOLDS:
    case D_TURB2LAM: case D_DESIGN:
      /*--- Already done in postprocessing of config ---*/
      break;
    default:
      break;
  }


  /*--- Initialize fan face pressure, fan face mach number, and mass flow rate ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inflow_MassFlow[iMarker]     = 0.0;
    Inflow_Mach[iMarker]         = Mach_Inf;
    Inflow_Pressure[iMarker]     = Pressure_Inf;
    Inflow_Area[iMarker]         = 0.0;
    
    Exhaust_MassFlow[iMarker]    = 0.0;
    Exhaust_Temperature[iMarker] = Temperature_Inf;
    Exhaust_Pressure[iMarker]    = Pressure_Inf;
    Exhaust_Area[iMarker]        = 0.0;
  }
  
  /*--- Initializate quantities for SlidingMesh Interface ---*/
  
  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++){
    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];

      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];

        SlidingStateNodes[iMarker][iPoint] = 0;
        for (iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }

    }
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CEulerVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/

  counter_local = 0;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    Density = node[iPoint]->GetSolution(0);

    Velocity2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);

    StaticEnergy= node[iPoint]->GetSolution(nDim+1)/Density - 0.5*Velocity2;

    FluidModel->SetTDState_rhoe(Density, StaticEnergy);
    Pressure= FluidModel->GetPressure();
    Temperature= FluidModel->GetTemperature();

    /*--- Use the values at the infinity ---*/

    if ((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
      Solution[0] = Density_Inf;
      for (iDim = 0; iDim < nDim; iDim++)
        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
      Solution[nDim+1] = Energy_Inf*Density_Inf;
      node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetSolution_Old(Solution);
      counter_local++;
    }

  }

  /*--- Warning message about non-physical points ---*/

  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  }

  /*--- Define solver parameters needed for execution of destructor ---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ) space_centered = true;
  else space_centered = false;

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;

  /*--- Perform the MPI communication of the solution ---*/

  Set_MPI_Solution(geometry, config);
  
}

CEulerSolver::~CEulerSolver(void) {

  unsigned short iVar, iMarker, iSpan;

  unsigned long iVertex;

  /*--- Array deallocation ---*/

  if (CD_Inv != NULL)         delete [] CD_Inv;
  if (CL_Inv != NULL)         delete [] CL_Inv;
  if (CSF_Inv != NULL)    delete [] CSF_Inv;
  if (CMx_Inv != NULL)           delete [] CMx_Inv;
  if (CMy_Inv != NULL)           delete [] CMy_Inv;
  if (CMz_Inv != NULL)           delete [] CMz_Inv;
  if (CFx_Inv != NULL)           delete [] CFx_Inv;
  if (CFy_Inv != NULL)           delete [] CFy_Inv;
  if (CFz_Inv != NULL)           delete [] CFz_Inv;
  if (CoPx_Inv != NULL)           delete [] CoPx_Inv;
  if (CoPy_Inv != NULL)           delete [] CoPy_Inv;
  if (CoPz_Inv != NULL)           delete [] CoPz_Inv;
  if (Surface_CL_Inv != NULL) delete[] Surface_CL_Inv;
  if (Surface_CD_Inv != NULL) delete[] Surface_CD_Inv;
  if (Surface_CSF_Inv != NULL) delete[] Surface_CSF_Inv;
  if (Surface_CEff_Inv != NULL) delete[] Surface_CEff_Inv;
  if (Surface_CFx_Inv != NULL)  delete [] Surface_CFx_Inv;
  if (Surface_CFy_Inv != NULL)  delete [] Surface_CFy_Inv;
  if (Surface_CFz_Inv != NULL)  delete [] Surface_CFz_Inv;
  if (Surface_CMx_Inv != NULL)  delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv != NULL)  delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv != NULL)  delete [] Surface_CMz_Inv;

  if (CD_Mnt != NULL)         delete [] CD_Mnt;
  if (CL_Mnt != NULL)         delete [] CL_Mnt;
  if (CSF_Mnt != NULL)    delete [] CSF_Mnt;
  if (CFx_Mnt != NULL)           delete [] CFx_Mnt;
  if (CFy_Mnt != NULL)           delete [] CFy_Mnt;
  if (CFz_Mnt != NULL)           delete [] CFz_Mnt;
  if (CMx_Mnt != NULL)           delete [] CMx_Mnt;
  if (CMy_Mnt != NULL)           delete [] CMy_Mnt;
  if (CMz_Mnt != NULL)           delete [] CMz_Mnt;
  if (CoPx_Mnt != NULL)           delete [] CoPx_Mnt;
  if (CoPy_Mnt != NULL)           delete [] CoPy_Mnt;
  if (CoPz_Mnt != NULL)           delete [] CoPz_Mnt;
  if (Surface_CL_Mnt != NULL) delete[] Surface_CL_Mnt;
  if (Surface_CD_Mnt != NULL) delete[] Surface_CD_Mnt;
  if (Surface_CSF_Mnt != NULL) delete[] Surface_CSF_Mnt;
  if (Surface_CEff_Mnt != NULL) delete[] Surface_CEff_Mnt;
  if (Surface_CFx_Mnt != NULL)  delete [] Surface_CFx_Mnt;
  if (Surface_CFy_Mnt != NULL)  delete [] Surface_CFy_Mnt;
  if (Surface_CFz_Mnt != NULL)  delete [] Surface_CFz_Mnt;
  if (Surface_CMx_Mnt != NULL)  delete [] Surface_CMx_Mnt;
  if (Surface_CMy_Mnt != NULL)  delete [] Surface_CMy_Mnt;
  if (Surface_CMz_Mnt != NULL)  delete [] Surface_CMz_Mnt;

  if (Surface_CL != NULL)    delete [] Surface_CL;
  if (Surface_CD != NULL)    delete [] Surface_CD;
  if (Surface_CSF != NULL) delete [] Surface_CSF;
  if (Surface_CEff != NULL) delete [] Surface_CEff;
  if (Surface_CFx != NULL)      delete [] Surface_CFx;
  if (Surface_CFy != NULL)      delete [] Surface_CFy;
  if (Surface_CFz != NULL)      delete [] Surface_CFz;
  if (Surface_CMx != NULL)      delete [] Surface_CMx;
  if (Surface_CMy != NULL)      delete [] Surface_CMy;
  if (Surface_CMz != NULL)      delete [] Surface_CMz;
  if (CEff_Inv != NULL)          delete [] CEff_Inv;
  if (CMerit_Inv != NULL)        delete [] CMerit_Inv;
  if (CT_Inv != NULL)            delete [] CT_Inv;
  if (CQ_Inv != NULL)            delete [] CQ_Inv;
  if (CEquivArea_Inv != NULL)    delete [] CEquivArea_Inv;
  if (CNearFieldOF_Inv != NULL)  delete [] CNearFieldOF_Inv;
  
  if (CEff_Mnt != NULL)          delete [] CEff_Mnt;
  if (CMerit_Mnt != NULL)        delete [] CMerit_Mnt;
  if (CT_Mnt != NULL)            delete [] CT_Mnt;
  if (CQ_Mnt != NULL)            delete [] CQ_Mnt;

  if (ForceInviscid != NULL)     delete [] ForceInviscid;
  if (MomentInviscid != NULL)    delete [] MomentInviscid;
  if (ForceMomentum != NULL)     delete [] ForceMomentum;
  if (MomentMomentum != NULL)    delete [] MomentMomentum;
  if (Inflow_MassFlow != NULL)  delete [] Inflow_MassFlow;
  if (Exhaust_MassFlow != NULL)  delete [] Exhaust_MassFlow;
  if (Exhaust_Area != NULL)      delete [] Exhaust_Area;
  if (Inflow_Pressure != NULL)  delete [] Inflow_Pressure;
  if (Inflow_Mach != NULL)      delete [] Inflow_Mach;
  if (Inflow_Area != NULL)      delete [] Inflow_Area;

  if (Exhaust_Pressure != NULL)  delete [] Exhaust_Pressure;
  if (Exhaust_Temperature != NULL)      delete [] Exhaust_Temperature;

  if (iPoint_UndLapl != NULL)       delete [] iPoint_UndLapl;
  if (jPoint_UndLapl != NULL)       delete [] jPoint_UndLapl;

  if (Primitive != NULL)        delete [] Primitive;
  if (Primitive_i != NULL)      delete [] Primitive_i;
  if (Primitive_j != NULL)      delete [] Primitive_j;

  if (Secondary != NULL)        delete [] Secondary;
  if (Secondary_i != NULL)      delete [] Secondary_i;
  if (Secondary_j != NULL)      delete [] Secondary_j;

  if (LowMach_Precontioner != NULL) {
    for (iVar = 0; iVar < nVar; iVar ++)
      delete [] LowMach_Precontioner[iVar];
    delete [] LowMach_Precontioner;
  }
  
  if (CPressure != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] CPressure[iMarker];
    delete [] CPressure;
  }
  
  if (CPressureTarget != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] CPressureTarget[iMarker];
    delete [] CPressureTarget;
  }
  
  if (CharacPrimVar != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] CharacPrimVar[iMarker][iVertex];
      delete [] CharacPrimVar[iMarker];
    }
    delete [] CharacPrimVar;
  }

  if (SlidingState != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if ( SlidingState[iMarker] != NULL ) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          if ( SlidingState[iMarker][iVertex] != NULL ){
            for (iVar = 0; iVar < nPrimVar+1; iVar++)
              delete [] SlidingState[iMarker][iVertex][iVar];
            delete [] SlidingState[iMarker][iVertex];
          }
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }
  
  if ( SlidingStateNodes != NULL ){
    for (iMarker = 0; iMarker < nMarker; iMarker++){
        if (SlidingStateNodes[iMarker] != NULL)
            delete [] SlidingStateNodes[iMarker];  
    }
    delete [] SlidingStateNodes;
  }
  
  if (DonorPrimVar != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] DonorPrimVar[iMarker][iVertex];
      delete [] DonorPrimVar[iMarker];
    }
    delete [] DonorPrimVar;
  }
  
  if (DonorGlobalIndex != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] DonorGlobalIndex[iMarker];
    delete [] DonorGlobalIndex;
  }

  if (ActDisk_DeltaP != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ActDisk_DeltaP[iMarker];
    delete [] ActDisk_DeltaP;
  }

  if (ActDisk_DeltaT != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ActDisk_DeltaT[iMarker];
    delete [] ActDisk_DeltaT;
  }

  if (Inlet_Ttotal != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Inlet_Ttotal[iMarker] != NULL)
        delete [] Inlet_Ttotal[iMarker];
    delete [] Inlet_Ttotal;
  }

  if (Inlet_Ptotal != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Inlet_Ptotal[iMarker] != NULL)
        delete [] Inlet_Ptotal[iMarker];
    delete [] Inlet_Ptotal;
  }

  if (Inlet_FlowDir != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Inlet_FlowDir[iMarker] != NULL) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          delete [] Inlet_FlowDir[iMarker][iVertex];
        delete [] Inlet_FlowDir[iMarker];
      }
    }
    delete [] Inlet_FlowDir;
  }

  if (nVertex != NULL)  delete [] nVertex;

  if (HeatFlux != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] HeatFlux[iMarker];
    }
    delete [] HeatFlux;
  }
  
  if (HeatFluxTarget != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] HeatFluxTarget[iMarker];
    }
    delete [] HeatFluxTarget;
  }
  
  if (YPlus != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] YPlus[iMarker];
    }
    delete [] YPlus;
  }

  if (Cauchy_Serie != NULL)  delete [] Cauchy_Serie;
  
  if (FluidModel != NULL) delete FluidModel;

  if(AverageVelocity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] AverageVelocity[iMarker][iSpan];
      delete  [] AverageVelocity[iMarker];
    }
    delete [] AverageVelocity;
  }

  if(AverageTurboVelocity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] AverageTurboVelocity[iMarker][iSpan];
      delete  [] AverageTurboVelocity[iMarker];
    }
    delete [] AverageTurboVelocity;
  }

  if(ExtAverageTurboVelocity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] ExtAverageTurboVelocity[iMarker][iSpan];
      delete  [] ExtAverageTurboVelocity[iMarker];
    }
    delete [] ExtAverageTurboVelocity;
  }


  if(AverageFlux !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] AverageFlux[iMarker][iSpan];
      delete  [] AverageFlux[iMarker];
    }
    delete [] AverageFlux;
  }

  if(SpanTotalFlux !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] SpanTotalFlux[iMarker][iSpan];
      delete  [] SpanTotalFlux[iMarker];
    }
    delete [] SpanTotalFlux;
  }

  if(AveragePressure !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] AveragePressure[iMarker];
    delete [] AveragePressure;
  }

  if(RadialEquilibriumPressure !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] RadialEquilibriumPressure[iMarker];
    delete [] RadialEquilibriumPressure;
  }

  if(ExtAveragePressure !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ExtAveragePressure[iMarker];
    delete [] ExtAveragePressure;
  }

  if(AverageDensity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] AverageDensity[iMarker];
    delete [] AverageDensity;
  }

  if(ExtAverageDensity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ExtAverageDensity[iMarker];
    delete [] ExtAverageDensity;
  }

  if(AverageKine !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] AverageKine[iMarker];
    delete [] AverageKine;
  }

  if(AverageOmega !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] AverageOmega[iMarker];
    delete [] AverageOmega;
  }

  if(AverageNu !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] AverageNu[iMarker];
    delete [] AverageNu;
  }

  if(ExtAverageKine !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ExtAverageKine[iMarker];
    delete [] ExtAverageKine;
  }

  if(ExtAverageOmega !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ExtAverageOmega[iMarker];
    delete [] ExtAverageOmega;
  }

  if(ExtAverageNu !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ExtAverageNu[iMarker];
    delete [] ExtAverageNu;
  }

  if(TurboVelocityIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
      for (iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++){
        delete [] TurboVelocityIn[iMarker][iSpan];
      }
      delete [] TurboVelocityIn[iMarker];
    }
    delete [] TurboVelocityIn;
  }

  if(TurboVelocityOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
      for (iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++){
        delete [] TurboVelocityOut[iMarker][iSpan];
      }
      delete [] TurboVelocityOut[iMarker];
    }
    delete [] TurboVelocityOut;
  }

  if(DensityIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] DensityIn[iMarker];
    delete [] DensityIn;
  }

  if(PressureIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] PressureIn[iMarker];
    delete [] PressureIn;
  }

  if(DensityOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] DensityOut[iMarker];
    delete [] DensityOut;
  }

  if(PressureOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] PressureOut[iMarker];
    delete [] PressureOut;
  }

  if(KineIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] KineIn[iMarker];
    delete [] KineIn;
  }

  if(OmegaIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] OmegaIn[iMarker];
    delete [] OmegaIn;
  }

  if(NuIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] NuIn[iMarker];
    delete [] NuIn;
  }

  if(KineOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] KineOut[iMarker];
    delete [] KineOut;
  }

  if(OmegaOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] OmegaOut[iMarker];
    delete [] OmegaOut;
  }

  if(NuOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] NuOut[iMarker];
    delete [] NuOut;
  }

  if(CkInflow !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan <nSpanWiseSections ; iSpan++)
        delete [] CkInflow[iMarker][iSpan];
      delete  [] CkInflow[iMarker];
    }
    delete [] CkInflow;
  }

  if(CkOutflow1 !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan <nSpanWiseSections ; iSpan++)
        delete [] CkOutflow1[iMarker][iSpan];
      delete  [] CkOutflow1[iMarker];
    }
    delete [] CkOutflow1;
  }

  if(CkOutflow2 !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan <nSpanWiseSections ; iSpan++)
        delete [] CkOutflow2[iMarker][iSpan];
      delete  [] CkOutflow2[iMarker];
    }
    delete [] CkOutflow2;
  }

}

void CEulerSolver::InitTurboContainers(CGeometry *geometry, CConfig *config){
  unsigned short iMarker, iSpan, iDim, iVar;
  unsigned short nSpanMax    = config->GetnSpanMaxAllZones();


  /*--- Initialize quantities for the average process for internal flow ---*/

  nSpanWiseSections = config->GetnSpanWiseSections();

  AverageVelocity                       = new su2double** [nMarker];
  AverageTurboVelocity                  = new su2double** [nMarker];
  OldAverageTurboVelocity               = new su2double** [nMarker];
  ExtAverageTurboVelocity               = new su2double** [nMarker];
  AverageFlux                           = new su2double** [nMarker];
  SpanTotalFlux                         = new su2double** [nMarker];
  AveragePressure                       = new su2double* [nMarker];
  OldAveragePressure                    = new su2double* [nMarker];
  RadialEquilibriumPressure             = new su2double* [nMarker];
  ExtAveragePressure                    = new su2double* [nMarker];
  AverageDensity                        = new su2double* [nMarker];
  OldAverageDensity                     = new su2double* [nMarker];
  ExtAverageDensity                     = new su2double* [nMarker];
  AverageNu                             = new su2double* [nMarker];
  AverageKine                           = new su2double* [nMarker];
  AverageOmega                          = new su2double* [nMarker];
  ExtAverageNu                          = new su2double* [nMarker];
  ExtAverageKine                        = new su2double* [nMarker];
  ExtAverageOmega                       = new su2double* [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    AverageVelocity[iMarker]                = new su2double* [nSpanWiseSections + 1];
    AverageTurboVelocity[iMarker]           = new su2double* [nSpanWiseSections + 1];
    OldAverageTurboVelocity[iMarker]        = new su2double* [nSpanWiseSections + 1];
    ExtAverageTurboVelocity[iMarker]        = new su2double* [nSpanWiseSections + 1];
    AverageFlux[iMarker]                    = new su2double* [nSpanWiseSections + 1];
    SpanTotalFlux[iMarker]                  = new su2double* [nSpanWiseSections + 1];
    AveragePressure[iMarker]                = new su2double [nSpanWiseSections + 1];
    OldAveragePressure[iMarker]             = new su2double [nSpanWiseSections + 1];
    RadialEquilibriumPressure[iMarker]      = new su2double [nSpanWiseSections + 1];
    ExtAveragePressure[iMarker]             = new su2double [nSpanWiseSections + 1];
    AverageDensity[iMarker]                 = new su2double [nSpanWiseSections + 1];
    OldAverageDensity[iMarker]              = new su2double [nSpanWiseSections + 1];
    ExtAverageDensity[iMarker]              = new su2double [nSpanWiseSections + 1];
    AverageNu[iMarker]                      = new su2double [nSpanWiseSections + 1];
    AverageKine[iMarker]                    = new su2double [nSpanWiseSections + 1];
    AverageOmega[iMarker]                   = new su2double [nSpanWiseSections + 1];
    ExtAverageNu[iMarker]                   = new su2double [nSpanWiseSections + 1];
    ExtAverageKine[iMarker]                 = new su2double [nSpanWiseSections + 1];
    ExtAverageOmega[iMarker]                = new su2double [nSpanWiseSections + 1];

    for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++){
      AverageVelocity[iMarker][iSpan]           = new su2double [nDim];
      AverageTurboVelocity[iMarker][iSpan]      = new su2double [nDim];
      OldAverageTurboVelocity[iMarker][iSpan]   = new su2double [nDim];
      ExtAverageTurboVelocity[iMarker][iSpan]   = new su2double [nDim];
      AverageFlux[iMarker][iSpan]               = new su2double [nVar];
      SpanTotalFlux[iMarker][iSpan]             = new su2double [nVar];
      AveragePressure[iMarker][iSpan]           = 0.0;
      OldAveragePressure[iMarker][iSpan]        = 0.0;
      RadialEquilibriumPressure[iMarker][iSpan] = 0.0;
      ExtAveragePressure[iMarker][iSpan]        = 0.0;
      AverageDensity[iMarker][iSpan]            = 0.0;
      OldAverageDensity[iMarker][iSpan]         = 0.0;
      ExtAverageDensity[iMarker][iSpan]         = 0.0;
      AverageNu[iMarker][iSpan]                 = 0.0;
      AverageKine[iMarker][iSpan]               = 0.0;
      AverageOmega[iMarker][iSpan]              = 0.0;
      ExtAverageNu[iMarker][iSpan]              = 0.0;
      ExtAverageKine[iMarker][iSpan]            = 0.0;
      ExtAverageOmega[iMarker][iSpan]           = 0.0;

      for (iDim = 0; iDim < nDim; iDim++) {
        AverageVelocity[iMarker][iSpan][iDim]            = 0.0;
        AverageTurboVelocity[iMarker][iSpan][iDim]       = 0.0;
        OldAverageTurboVelocity[iMarker][iSpan][iDim]    = 0.0;
        ExtAverageTurboVelocity[iMarker][iSpan][iDim]    = 0.0;
      }

      for (iVar = 0; iVar < nVar; iVar++) {
        AverageFlux[iMarker][iSpan][iVar]       = 0.0;
        SpanTotalFlux[iMarker][iSpan][iVar]     = 0.0;
      }
    }
  }


  /*--- Initialize primitive quantities for turboperformace ---*/

  nMarkerTurboPerf = config->GetnMarker_TurboPerformance();

  DensityIn                     = new su2double*[nMarkerTurboPerf];
  PressureIn                    = new su2double*[nMarkerTurboPerf];
  TurboVelocityIn               = new su2double**[nMarkerTurboPerf];
  DensityOut                    = new su2double*[nMarkerTurboPerf];
  PressureOut                   = new su2double*[nMarkerTurboPerf];
  TurboVelocityOut              = new su2double**[nMarkerTurboPerf];
  KineIn                        = new su2double*[nMarkerTurboPerf];
  OmegaIn                       = new su2double*[nMarkerTurboPerf];
  NuIn                          = new su2double*[nMarkerTurboPerf];
  KineOut                       = new su2double*[nMarkerTurboPerf];
  OmegaOut                      = new su2double*[nMarkerTurboPerf];
  NuOut                         = new su2double*[nMarkerTurboPerf];

  for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
    DensityIn[iMarker]          = new su2double [nSpanMax + 1];
    PressureIn[iMarker]         = new su2double [nSpanMax + 1];
    TurboVelocityIn[iMarker]    = new su2double*[nSpanMax + 1];
    DensityOut[iMarker]         = new su2double [nSpanMax + 1];
    PressureOut[iMarker]        = new su2double [nSpanMax + 1];
    TurboVelocityOut[iMarker]   = new su2double*[nSpanMax + 1];
    KineIn[iMarker]             = new su2double [nSpanMax + 1];
    OmegaIn[iMarker]            = new su2double [nSpanMax + 1];
    NuIn[iMarker]               = new su2double [nSpanMax + 1];
    KineOut[iMarker]            = new su2double [nSpanMax + 1];
    OmegaOut[iMarker]           = new su2double [nSpanMax + 1];
    NuOut[iMarker]              = new su2double [nSpanMax + 1];



    for (iSpan = 0; iSpan < nSpanMax + 1; iSpan++){
      DensityIn[iMarker][iSpan]          = 0.0;
      PressureIn[iMarker][iSpan]         = 0.0;
      TurboVelocityIn[iMarker][iSpan]    = new su2double[nDim];
      DensityOut[iMarker][iSpan]         = 0.0;
      PressureOut[iMarker][iSpan]        = 0.0;
      TurboVelocityOut [iMarker][iSpan]  = new su2double[nDim];
      KineIn[iMarker][iSpan]             = 0.0;
      OmegaIn[iMarker][iSpan]            = 0.0;
      NuIn[iMarker][iSpan]               = 0.0;
      KineOut[iMarker][iSpan]            = 0.0;
      OmegaOut[iMarker][iSpan]           = 0.0;
      NuOut[iMarker][iSpan]              = 0.0;

      for (iDim = 0; iDim < nDim; iDim++){
        TurboVelocityIn  [iMarker][iSpan][iDim]   = 0.0;
        TurboVelocityOut [iMarker][iSpan][iDim]   = 0.0;
      }
    }
  }


  /*--- Initialize quantities for NR BC ---*/

  if(config->GetBoolGiles()){

    CkInflow= new complex<su2double>** [nMarker];
    CkOutflow1= new complex<su2double>** [nMarker];
    CkOutflow2= new complex<su2double>** [nMarker];

    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      CkInflow[iMarker]= new complex<su2double>*[nSpanWiseSections];
      CkOutflow1[iMarker]= new complex<su2double>*[nSpanWiseSections];
      CkOutflow2[iMarker]= new complex<su2double>*[nSpanWiseSections];

      for(iSpan = 0; iSpan < nSpanWiseSections; iSpan++) {
        CkInflow[iMarker][iSpan]= new complex<su2double>[2*geometry->GetnFreqSpanMax(INFLOW)+1];
        CkOutflow1[iMarker][iSpan]= new complex<su2double>[2*geometry->GetnFreqSpanMax(OUTFLOW)+1];
        CkOutflow2[iMarker][iSpan]= new complex<su2double>[2*geometry->GetnFreqSpanMax(OUTFLOW)+1];

        for (iVar = 0; iVar <2*geometry->GetnFreqSpanMax(INFLOW)+1; iVar++) {
          CkInflow[iMarker][iSpan][iVar]= complex<su2double>(0.0,0.0);
        }

        for (iVar = 0; iVar <2*geometry->GetnFreqSpanMax(OUTFLOW)+1; iVar++) {
          CkOutflow1[iMarker][iSpan][iVar]= complex<su2double>(0.0,0.0);
          CkOutflow2[iMarker][iSpan][iVar]= complex<su2double>(0.0,0.0);
        }
      }
    }
  }


}

void CEulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CEulerSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
}

void CEulerSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Undivided_Laplacian = new su2double [nBufferR_Vector];
      Buffer_Send_Undivided_Laplacian = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Undivided_Laplacian[iVar*nVertexS+iVertex] = node[iPoint]->GetUndivided_Laplacian(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex] = Buffer_Send_Undivided_Laplacian[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Undivided_Laplacian;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetUndivided_Laplacian(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Undivided_Laplacian;
      
    }
    
  }
  
}

void CEulerSolver::Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR, *Buffer_Receive_Neighbor = NULL, *Buffer_Send_Neighbor = NULL;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new su2double [nBufferR_Vector];
      Buffer_Send_Lambda = new su2double[nBufferS_Vector];
      Buffer_Receive_Neighbor = new unsigned short [nBufferR_Vector];
      Buffer_Send_Neighbor = new unsigned short[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetLambda();
        Buffer_Send_Neighbor[iVertex] = geometry->node[iPoint]->GetnPoint();
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      SU2_MPI::Sendrecv(Buffer_Send_Neighbor, nBufferS_Vector, MPI_UNSIGNED_SHORT, send_to, 1,
                        Buffer_Receive_Neighbor, nBufferR_Vector, MPI_UNSIGNED_SHORT, receive_from, 1, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
        Buffer_Receive_Neighbor[iVertex] = Buffer_Send_Neighbor[iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      delete [] Buffer_Send_Neighbor;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetLambda(Buffer_Receive_Lambda[iVertex]);
        geometry->node[iPoint]->SetnNeighbor(Buffer_Receive_Neighbor[iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      delete [] Buffer_Receive_Neighbor;
      
    }
    
  }
}

void CEulerSolver::Set_MPI_Sensor(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new su2double [nBufferR_Vector];
      Buffer_Send_Lambda = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetSensor();
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetSensor(Buffer_Receive_Lambda[iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      
    }
    
  }
}

void CEulerSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CEulerSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nVar];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetLimiter(iVar, Limiter[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
  }
  
  delete [] Limiter;
  
}

void CEulerSolver::Set_MPI_Primitive_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad*nDim;        nBufferR_Vector = nVertexR*nPrimVarGrad*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient_Primitive(iVar, iDim);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient_Primitive(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
  }
  
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CEulerSolver::Set_MPI_Primitive_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nPrimVarGrad];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad;        nBufferR_Vector = nVertexR*nPrimVarGrad;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter_Primitive(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          node[iPoint]->SetLimiter_Primitive(iVar, Limiter[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
  }
  
  delete [] Limiter;
  
}

void CEulerSolver::Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0;
  long iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker;
  long nDomain = 0, iDomain, jDomain;
  //bool ActDisk_Perimeter;
  bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS));
  
  unsigned short nPrimVar_ = nPrimVar;
  if (rans) nPrimVar_ += 2; // Add two extra variables for the turbulence.
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status;
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double *Buffer_Send_PrimVar = NULL;
  long      *Buffer_Send_Data    = NULL;
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  su2double *iPrimVar = new su2double [nPrimVar_];
  
  unsigned long Buffer_Size_PrimVar = 0;
  unsigned long Buffer_Size_Data    = 0;
  
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double *Buffer_Receive_PrimVar = NULL;
  long      *Buffer_Receive_Data    = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          //ActDisk_Perimeter = geometry->vertex[iMarker][iVertex]->GetActDisk_Perimeter();
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
//          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain()) && (!ActDisk_Perimeter)) {
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            Buffer_Send_nPointTotal++;
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar += nPointTotal_s[iDomain]*(nPrimVar_);
    Buffer_Size_Data += nPointTotal_s[iDomain]*(3);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar = new su2double[Buffer_Size_PrimVar];
  Buffer_Send_Data    = new long[Buffer_Size_Data];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          //ActDisk_Perimeter = geometry->vertex[iMarker][iVertex]->GetActDisk_Perimeter();
          
//          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain()) && (!ActDisk_Perimeter)) {
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              Buffer_Send_PrimVar[(nPrimVar_)*(PointTotal_Counter+iPointTotal)+iVar] = node[iPoint]->GetPrimitive(iVar);
            }
            if (rans) {
              Buffer_Send_PrimVar[(nPrimVar_)*(PointTotal_Counter+iPointTotal)+nPrimVar] = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
              Buffer_Send_PrimVar[(nPrimVar_)*(PointTotal_Counter+iPointTotal)+(nPrimVar+1)] = 0.0;
            }
            
            iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
            jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
            jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
            
            Buffer_Send_Data[(3)*(PointTotal_Counter+iPointTotal)+(0)]  = iGlobalIndex;
            Buffer_Send_Data[(3)*(PointTotal_Counter+iPointTotal)+(1)] = jVertex;
            Buffer_Send_Data[(3)*(PointTotal_Counter+iPointTotal)+(2)]  = jMarker;
            
            iPointTotal++;
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar_)],
                     nPointTotal_s[iDomain]*(nPrimVar_), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
      SU2_MPI::Bsend(&Buffer_Send_Data[PointTotal_Counter*(3)],
                     nPointTotal_s[iDomain]*(3), MPI_LONG, iDomain,
                     iDomain+nDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(nPrimVar_)];
      Buffer_Receive_Data               = new long[nPointTotal_s[iDomain]*(3)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(nPrimVar_); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar_)+iter];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(3); iter++)
        Buffer_Receive_Data[iter] = Buffer_Send_Data[PointTotal_Counter*(3)+iter];
      
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        for (iVar = 0; iVar < nPrimVar_; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar_)+iVar];
        
        iGlobal       =  Buffer_Receive_Data[iPoint*(3)+(0)];
        iVertex      =  Buffer_Receive_Data[iPoint*(3)+(1)];
        iMarker      = Buffer_Receive_Data[iPoint*(3)+(2)];
        
        for (iVar = 0; iVar < nPrimVar_; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      delete [] Buffer_Receive_Data;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(nPrimVar_)];
      Buffer_Receive_Data               = new long [nPointTotal_r[iDomain]*(3)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(nPrimVar_) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status);
      
      SU2_MPI::Recv(Buffer_Receive_Data, nPointTotal_r[iDomain]*(3) , MPI_LONG,
                    iDomain, rank+nDomain, MPI_COMM_WORLD, &status);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = Buffer_Receive_Data[iPoint*(3)+(0)];
        iVertex      = Buffer_Receive_Data[iPoint*(3)+(1)];
        iMarker      = Buffer_Receive_Data[iPoint*(3)+(2)];
        
        for (iVar = 0; iVar < nPrimVar_; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar_)+iVar];
        
        for (iVar = 0; iVar < nPrimVar_; iVar++) {
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        }
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      delete [] Buffer_Receive_Data;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  delete[] Buffer_Send_Data;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
}

void CEulerSolver::Set_MPI_Nearfield(CGeometry *geometry, CConfig *config) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0;
  long iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker;
  long nDomain = 0, iDomain, jDomain;
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  

#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  su2double        *iPrimVar          = new su2double [nPrimVar];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            Buffer_Send_nPointTotal++;
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(nPrimVar+3);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
            jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
            jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+iVar] = node[iPoint]->GetPrimitive(iVar);
            }
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+0)]  = su2double(iGlobalIndex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+1)] = su2double(jVertex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+2)]  = su2double(jMarker);
            
            iPointTotal++;
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)],
                     nPointTotal_s[iDomain]*(nPrimVar+3), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(nPrimVar+3)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(nPrimVar+3); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal       =  SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      =  SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(nPrimVar+3)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(nPrimVar+3) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar,  iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
}

void CEulerSolver::Set_MPI_Interface(CGeometry *geometry, CConfig *config) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker;
  long nDomain = 0, iDomain, jDomain;
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  

#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [nPrimVar];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            Buffer_Send_nPointTotal++;
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(nPrimVar+3);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            
            iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
            jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
            jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+iVar] = node[iPoint]->GetPrimitive(iVar);
            }
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+0)]  = su2double(iGlobalIndex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+1)] = su2double(jVertex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+2)]  = su2double(jMarker);
            
            iPointTotal++;
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)],
                     nPointTotal_s[iDomain]*(nPrimVar+3), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(nPrimVar+3)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(nPrimVar+3); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal       =  SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(nPrimVar+3)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(nPrimVar+3) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
}

void CEulerSolver::SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh) {
  
  su2double Temperature_FreeStream = 0.0, Mach2Vel_FreeStream = 0.0, ModVel_FreeStream = 0.0,
  Energy_FreeStream = 0.0, ModVel_FreeStreamND = 0.0, Velocity_Reynolds = 0.0,
  Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Tke_FreeStream = 0.0,
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
  Temperature_Ref = 0.0, Time_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Energy_Ref= 0.0,
  Froude = 0.0, Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0,
  Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0, TgammaR = 0.0;

  unsigned short iDim;

  /*--- Local variables ---*/
  
  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach             = config->GetMach();
  su2double Reynolds         = config->GetReynolds();
  bool unsteady           = (config->GetUnsteady_Simulation() != NO);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool gravity            = config->GetGravityForce();
  bool turbulent          = (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS);
  bool tkeNeeded          = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == TEMPERATURE_FS);
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);
  bool aeroelastic        = config->GetAeroelastic_Simulation();
  
  /*--- Set temperature via the flutter speed index ---*/
  if (aeroelastic) {
    su2double vf             = config->GetAeroelastic_Flutter_Speed_Index();
    su2double w_alpha        = config->GetAeroelastic_Frequency_Pitch();
    su2double b              = config->GetLength_Reynolds()/2.0; // airfoil semichord, Reynolds length is by defaul 1.0
    su2double mu             = config->GetAeroelastic_Airfoil_Mass_Ratio();
    // The temperature times gamma times the gas constant. Depending on the FluidModel temp is calculated below.
    TgammaR = ((vf*vf)*(b*b)*(w_alpha*w_alpha)*mu) / (Mach*Mach);
  }
  
  /*--- Compressible non dimensionalization ---*/

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream  = config->GetTemperature_FreeStream();

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:

      if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
      else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(1716.49);

      FluidModel = new CIdealGas(1.4, config->GetGas_Constant());
      if (free_stream_temp) {
        if (aeroelastic) {
          Temperature_FreeStream = TgammaR / (config->GetGas_Constant()*1.4);
          config->SetTemperature_FreeStream(Temperature_FreeStream);
        }
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case IDEAL_GAS:

      FluidModel = new CIdealGas(Gamma, config->GetGas_Constant());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case VW_GAS:

      FluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(),
                                       config->GetPressure_Critical(), config->GetTemperature_Critical());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case PR_GAS:

      FluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(),
                                     config->GetTemperature_Critical(), config->GetAcentric_Factor());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

  }

  Mach2Vel_FreeStream = FluidModel->GetSoundSpeed();

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  if (nDim == 2) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
  }
  if (nDim == 3) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
  }

  /*--- Compute the modulus of the free stream velocity ---*/

  ModVel_FreeStream = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  /*--- Viscous initialization ---*/

  if (viscous) {

    /*--- The dimensional viscosity is needed to determine the free-stream conditions.
          To accomplish this, simply set the non-dimensional coefficients to the
          dimensional ones. This will be overruled later.---*/
    config->SetMu_RefND(config->GetMu_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref());
    config->SetMu_SND(config->GetMu_S());

    config->SetMu_ConstantND(config->GetMu_Constant());

    /*--- Reynolds based initialization ---*/

    if (reynolds_init) {

      /*--- First, check if there is mesh motion. If yes, use the Mach
         number relative to the body to initialize the flow. ---*/

      if (grid_movement) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
      else Velocity_Reynolds = ModVel_FreeStream;

      /*--- For viscous flows, pressure will be computed from a density
            that is found from the Reynolds number. The viscosity is computed
            from the dimensional version of Sutherland's law or the constant
            viscosity, depending on the input option.---*/

      FluidModel->SetLaminarViscosityModel(config);

      Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);

      Density_FreeStream = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*config->GetLength_Reynolds());
      config->SetDensity_FreeStream(Density_FreeStream);
      FluidModel->SetTDState_rhoT(Density_FreeStream, Temperature_FreeStream);
      Pressure_FreeStream = FluidModel->GetPressure();
      config->SetPressure_FreeStream(Pressure_FreeStream);
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Thermodynamics quantities based initialization ---*/

    else {

      FluidModel->SetLaminarViscosityModel(config);
      Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Turbulence kinetic energy ---*/

    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());

  }
  else {

    /*--- For inviscid flow, energy is calculated from the specified
       FreeStream quantities using the proper gas law. ---*/

    Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

  }

  /*-- Compute the freestream energy. ---*/

  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Compute non dimensional quantities. By definition,
     Lref is one because we have converted the grid to meters. ---*/

  if (config->GetRef_NonDim() == DIMENSIONAL) {
    Pressure_Ref      = 1.0;
    Density_Ref       = 1.0;
    Temperature_Ref   = 1.0;
  }
  else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
    Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
    Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
    Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
    Pressure_Ref      = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  config->SetPressure_Ref(Pressure_Ref);
  config->SetDensity_Ref(Density_Ref);
  config->SetTemperature_Ref(Temperature_Ref);

  Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
  Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
  Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
  Force_Ref         = config->GetDensity_Ref()*Velocity_Ref*Velocity_Ref*Length_Ref*Length_Ref; config->SetForce_Ref(Force_Ref);
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref  = Viscosity_Ref*Gas_Constant_Ref;                              config->SetConductivity_Ref(Conductivity_Ref);
  Froude            = ModVel_FreeStream/sqrt(STANDART_GRAVITY*Length_Ref);         config->SetFroude(Froude);
  
  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
  
  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();  config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();    config->SetDensity_FreeStreamND(Density_FreeStreamND);
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }
  
  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);
  
  Gas_ConstantND = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);
  
  
  ModVel_FreeStreamND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND); config->SetModVel_FreeStreamND(ModVel_FreeStreamND);
  
  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;   config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);
  
  Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStream(Tke_FreeStream);
  
  Tke_FreeStreamND  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStreamND(Tke_FreeStreamND);
  
  Omega_FreeStream = Density_FreeStream*Tke_FreeStream/(Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStream(Omega_FreeStream);
  
  Omega_FreeStreamND = Density_FreeStreamND*Tke_FreeStreamND/(Viscosity_FreeStreamND*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStreamND(Omega_FreeStreamND);
  
  /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/
 
  /*--- Delete the original (dimensional) FluidModel object before replacing. ---*/
  
  delete FluidModel;
 
  switch (config->GetKind_FluidModel()) {
      
    case STANDARD_AIR:
      FluidModel = new CIdealGas(1.4, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case IDEAL_GAS:
      FluidModel = new CIdealGas(Gamma, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case VW_GAS:
      FluidModel = new CVanDerWaalsGas(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                       config->GetTemperature_Critical()/config->GetTemperature_Ref());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case PR_GAS:
      FluidModel = new CPengRobinson(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                     config->GetTemperature_Critical()/config->GetTemperature_Ref(), config->GetAcentric_Factor());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
  }
  
  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
  
  if (viscous) {
    
    /*--- Constant viscosity model ---*/
    config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);
    
    /*--- Sutherland's model ---*/
    
    config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_S()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()/config->GetTemperature_Ref());
    
    /* constant thermal conductivity model */
    config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);
    
    FluidModel->SetLaminarViscosityModel(config);
    FluidModel->SetThermalConductivityModel(config);
    
  }
  
  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);
  
  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);
  
  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);
  
  /*--- Write output to the console if this is the master node and first domain ---*/
  
  if ((rank == MASTER_NODE) && (iMesh == MESH_0)) {
    
    cout.precision(6);
    
    if (viscous) {
      cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
      cout << "based on the free-stream temperature and a density computed" << endl;
      cout << "from the Reynolds number." << endl;
    } else {
      cout << "Inviscid flow: Computing density based on free-stream" << endl;
      cout << "temperature and pressure using the ideal gas law." << endl;
    }
    
    if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
    else cout << "Force coefficients computed using free-stream values." << endl;
    
    cout <<"-- Input conditions:"<< endl;
    
    switch (config->GetKind_FluidModel()) {

      case STANDARD_AIR:
        cout << "Fluid Model: STANDARD_AIR "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant();
        if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        break;

      case IDEAL_GAS:
        cout << "Fluid Model: IDEAL_GAS "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        break;

      case VW_GAS:
        cout << "Fluid Model: Van der Waals "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
        cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
        cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
        cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
        break;

      case PR_GAS:
        cout << "Fluid Model: Peng-Robinson "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
        cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
        cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
        cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
        break;

    }
    if (viscous) {
      switch (config->GetKind_ViscosityModel()) {

        case CONSTANT_VISCOSITY:
          cout << "Viscosity Model: CONSTANT_VISCOSITY  "<< endl;
          cout << "Laminar Viscosity: " << config->GetMu_Constant();
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          break;

        case SUTHERLAND:
          cout << "Viscosity Model: SUTHERLAND "<< endl;
          cout << "Ref. Laminar Viscosity: " << config->GetMu_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Ref. Temperature: " << config->GetMu_Temperature_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Sutherland Constant: "<< config->GetMu_S();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          cout << "Ref. Temperature (non-dim): " << config->GetMu_Temperature_RefND()<< endl;
          cout << "Sutherland constant (non-dim): "<< config->GetMu_SND()<< endl;
          break;

      }
      switch (config->GetKind_ConductivityModel()) {

        case CONSTANT_PRANDTL:
          cout << "Conductivity Model: CONSTANT_PRANDTL  "<< endl;
          cout << "Prandtl: " << config->GetPrandtl_Lam()<< endl;
          break;

        case CONSTANT_CONDUCTIVITY:
          cout << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< endl;
          cout << "Molecular Conductivity: " << config->GetKt_Constant()<< " W/m^2.K." << endl;
          cout << "Molecular Conductivity (non-dim): " << config->GetKt_ConstantND()<< endl;
          break;

      }
    }

    
    cout << "Free-stream static pressure: " << config->GetPressure_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Free-stream total pressure: " << config->GetPressure_FreeStream() * pow( 1.0+Mach*Mach*0.5*(Gamma-1.0), Gamma/(Gamma-1.0) );
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Free-stream temperature: " << config->GetTemperature_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
    
    cout << "Free-stream density: " << config->GetDensity_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;
    
    if (nDim == 2) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ")";
    }
    if (nDim == 3) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ", " << config->GetVelocity_FreeStream()[2] << ")";
    }
    if (config->GetSystemMeasurements() == SI) cout << " m/s. ";
    else if (config->GetSystemMeasurements() == US) cout << " ft/s. ";
    
    cout << "Magnitude: "   << config->GetModVel_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m/s (" << config->GetModVel_FreeStream()*1.94384 << " KTS)." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s (" << config->GetModVel_FreeStream()*0.592484 << " KTS)." << endl;
    
    cout << "Free-stream total energy per unit mass: " << config->GetEnergy_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
    
    if (viscous) {
      cout << "Free-stream viscosity: " << config->GetViscosity_FreeStream();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      if (turbulent) {
        cout << "Free-stream turb. kinetic energy per unit mass: " << config->GetTke_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
        cout << "Free-stream specific dissipation: " << config->GetOmega_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " 1/s." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " 1/s." << endl;
      }
    }
    
    if (unsteady) { cout << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime() << " s." << endl; }
    
    /*--- Print out reference values. ---*/
    
    cout <<"-- Reference values:"<< endl;
    
    cout << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;

    cout << "Reference pressure: " << config->GetPressure_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Reference temperature: " << config->GetTemperature_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
    
    cout << "Reference density: " << config->GetDensity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;
    
    cout << "Reference velocity: " << config->GetVelocity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;
    
    cout << "Reference energy per unit mass: " << config->GetEnergy_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
    
    if (viscous) {
      cout << "Reference viscosity: " << config->GetViscosity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      cout << "Reference conductivity: " << config->GetConductivity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " W/m^2.K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf/ft.s.R." << endl;
    }
    
    
    if (unsteady) cout << "Reference time: " << config->GetTime_Ref() <<" s." << endl;
    
    /*--- Print out resulting non-dim values here. ---*/
    
    cout << "-- Resulting non-dimensional state:" << endl;
    cout << "Mach number (non-dim): " << config->GetMach() << endl;
    if (viscous) {
      cout << "Reynolds number (non-dim): " << config->GetReynolds() <<". Re length: " << config->GetLength_Reynolds();
      if (config->GetSystemMeasurements() == SI) cout << " m." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " ft." << endl;
    }
    if (gravity) {
      cout << "Froude number (non-dim): " << Froude << endl;
      cout << "Lenght of the baseline wave (non-dim): " << 2.0*PI_NUMBER*Froude*Froude << endl;
    }
    
    cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << endl;
    cout << "Free-stream temperature (non-dim): " << config->GetTemperature_FreeStreamND() << endl;

    cout << "Free-stream pressure (non-dim): " << config->GetPressure_FreeStreamND() << endl;
    
    cout << "Free-stream density (non-dim): " << config->GetDensity_FreeStreamND() << endl;
    
    if (nDim == 2) {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
    }
    cout << "Magnitude: "    << config->GetModVel_FreeStreamND() << endl;
    
    cout << "Free-stream total energy per unit mass (non-dim): " << config->GetEnergy_FreeStreamND() << endl;
    
    if (viscous) {
      cout << "Free-stream viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << endl;
      if (turbulent) {
        cout << "Free-stream turb. kinetic energy (non-dim): " << config->GetTke_FreeStreamND() << endl;
        cout << "Free-stream specific dissipation (non-dim): " << config->GetOmega_FreeStreamND() << endl;
      }
    }
    
    if (unsteady) {
      cout << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << endl;
      cout << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << endl;
    }
    
    cout << endl;
    
  }
  
}

void CEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {

  unsigned long iPoint;
  unsigned short iMesh, iDim;
  su2double X0[3] = {0.0,0.0,0.0}, X1[3] = {0.0,0.0,0.0}, X2[3] = {0.0,0.0,0.0},
  X1_X0[3] = {0.0,0.0,0.0}, X2_X0[3] = {0.0,0.0,0.0}, X2_X1[3] = {0.0,0.0,0.0},
  CP[3] = {0.0,0.0,0.0}, Distance, DotCheck, Radius;
  
  unsigned short nDim = geometry[MESH_0]->GetnDim();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool rans = ((config->GetKind_Solver() == RANS) ||
               (config->GetKind_Solver() == ADJ_RANS) ||
               (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool SubsonicEngine = config->GetSubsonicEngine();

  /*--- Set subsonic initial condition for engine intakes ---*/
  
  if (SubsonicEngine) {
    
    /*--- Set initial boundary condition at iteration 0 ---*/
    
    if ((ExtIter == 0) && (!restart)) {
      
      su2double Velocity_Cyl[3] = {0.0, 0.0, 0.0}, Velocity_CylND[3] = {0.0, 0.0, 0.0}, Viscosity_Cyl,
      Density_Cyl, Density_CylND, Pressure_CylND, ModVel_Cyl, ModVel_CylND, Energy_CylND,
      T_ref = 0.0, S = 0.0, Mu_ref = 0.0, *Coord, *SubsonicEngine_Cyl, *SubsonicEngine_Values;
      
      SubsonicEngine_Values = config->GetSubsonicEngine_Values();
      su2double Mach_Cyl        = SubsonicEngine_Values[0];
      su2double Alpha_Cyl       = SubsonicEngine_Values[1];
      su2double Beta_Cyl        = SubsonicEngine_Values[2];
      su2double Pressure_Cyl    = SubsonicEngine_Values[3];
      su2double Temperature_Cyl = SubsonicEngine_Values[4];
      
      su2double Alpha = Alpha_Cyl*PI_NUMBER/180.0;
      su2double Beta  = Beta_Cyl*PI_NUMBER/180.0;
      
      su2double Gamma_Minus_One = Gamma - 1.0;
      su2double Gas_Constant = config->GetGas_Constant();
      
      su2double Mach2Vel_Cyl = sqrt(Gamma*Gas_Constant*Temperature_Cyl);
      
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          
          Velocity_Cyl[0] = cos(Alpha)*cos(Beta)*Mach_Cyl*Mach2Vel_Cyl;
          Velocity_Cyl[1] = sin(Beta)*Mach_Cyl*Mach2Vel_Cyl;
          Velocity_Cyl[2] = sin(Alpha)*cos(Beta)*Mach_Cyl*Mach2Vel_Cyl;
          
          ModVel_Cyl = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            ModVel_Cyl += Velocity_Cyl[iDim]*Velocity_Cyl[iDim];
          }
          ModVel_Cyl = sqrt(ModVel_Cyl);
          
          if (config->GetViscous()) {
            if (config->GetSystemMeasurements() == SI) { T_ref = 273.15; S = 110.4; Mu_ref = 1.716E-5; }
            if (config->GetSystemMeasurements() == US) {
              T_ref = (273.15 - 273.15) * 1.8 + 491.67;
              S = (110.4 - 273.15) * 1.8 + 491.67;
              Mu_ref = 1.716E-5/47.88025898;
            }
            Viscosity_Cyl = Mu_ref*(pow(Temperature_Cyl/T_ref, 1.5) * (T_ref+S)/(Temperature_Cyl+S));
            Density_Cyl   = config->GetReynolds()*Viscosity_Cyl/(ModVel_Cyl*config->GetLength_Reynolds());
            Pressure_Cyl  = Density_Cyl*Gas_Constant*Temperature_Cyl;
          }
          else {
            Density_Cyl = Pressure_Cyl/(Gas_Constant*Temperature_Cyl);
          }
          
          Density_CylND  = Density_Cyl/config->GetDensity_Ref();
          Pressure_CylND = Pressure_Cyl/config->GetPressure_Ref();
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_CylND[iDim] = Velocity_Cyl[iDim]/config->GetVelocity_Ref();
          }
          
          ModVel_CylND = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            ModVel_CylND += Velocity_CylND[iDim]*Velocity_CylND[iDim];
          }
          ModVel_CylND = sqrt(ModVel_CylND);
          
          Energy_CylND = Pressure_CylND/(Density_CylND*Gamma_Minus_One)+0.5*ModVel_CylND*ModVel_CylND;
          
          Coord = geometry[iMesh]->node[iPoint]->GetCoord();
          
          SubsonicEngine_Cyl = config->GetSubsonicEngine_Cyl();
          
          X0[0] = Coord[0];                             X0[1] = Coord[1];                           X0[2] = Coord[2];
          X1[0] = SubsonicEngine_Cyl[0]; X1[1] = SubsonicEngine_Cyl[1]; X1[2] = SubsonicEngine_Cyl[2];
          X2[0] = SubsonicEngine_Cyl[3]; X2[1] = SubsonicEngine_Cyl[4]; X2[2] = SubsonicEngine_Cyl[5];
          Radius = SubsonicEngine_Cyl[6];
          
          for (iDim = 0; iDim < nDim; iDim++) {
            X2_X1[iDim]= X1[iDim] - X2[iDim];
            X1_X0[iDim]= X0[iDim] - X1[iDim];
            X2_X0[iDim]= X0[iDim] - X2[iDim];
          }
          
          CP[0] = (X2_X1[1]*X1_X0[2] - X2_X1[2]*X1_X0[1]);
          CP[1] = (X2_X1[2]*X1_X0[0] - X2_X1[0]*X1_X0[2]);
          CP[2] = (X2_X1[0]*X1_X0[1] - X2_X1[1]*X1_X0[0]);
          
          Distance = sqrt((CP[0]*CP[0]+CP[1]*CP[1]+CP[2]*CP[2])/(X2_X1[0]*X2_X1[0]+X2_X1[1]*X2_X1[1]+X2_X1[2]*X2_X1[2]));
          
          DotCheck = -(X1_X0[0]*X2_X1[0]+X1_X0[1]*X2_X1[1]+X1_X0[2]*X2_X1[2]);
          if (DotCheck < 0.0) Distance = sqrt(X1_X0[0]*X1_X0[0]+X1_X0[1]*X1_X0[1]+X1_X0[2]*X1_X0[2]);
          
          DotCheck = (X2_X0[0]*X2_X1[0]+X2_X0[1]*X2_X1[1]+X2_X0[2]*X2_X1[2]);
          if (DotCheck < 0.0) Distance = sqrt(X2_X0[0]*X2_X0[0]+X2_X0[1]*X2_X0[1]+X2_X0[2]*X2_X0[2]);
          
          if (Distance < Radius) {
            
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(0, Density_CylND);
            for (iDim = 0; iDim < nDim; iDim++)
              solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(iDim+1, Density_CylND*Velocity_CylND[iDim]);
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(nVar-1, Density_CylND*Energy_CylND);
            
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(0, Density_CylND);
            for (iDim = 0; iDim < nDim; iDim++)
              solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(iDim+1, Density_CylND*Velocity_CylND[iDim]);
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(nVar-1, Density_CylND*Energy_CylND);
            
          }
          
        }
        
        /*--- Set the MPI communication ---*/
        
        solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
        solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution_Old(geometry[iMesh], config);
        
      }
      
    }
    
  }

  /*--- Make sure that the solution is well initialized for unsteady
   calculations with dual time-stepping (load additional restarts for 2nd-order). ---*/
  
  if (dual_time && (ExtIter == 0 || (restart && (long)ExtIter == config->GetUnst_RestartIter()))) {
    
    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/
    
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
        if (rans) {
          solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n1();
        }
      }
    }
    
    if ((restart && (long)ExtIter == config->GetUnst_RestartIter()) &&
        (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
      
      /*--- Load an additional restart file for a 2nd-order restart ---*/
      
      solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1), true);
      
      /*--- Load an additional restart file for the turbulence model ---*/
      if (rans)
        solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1), false);
      
      /*--- Push back this new solution to time level N. ---*/
      
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          if (rans) {
            solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          }
        }
      }
    }
  }

}

void CEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long ErrorCounter = 0;
  
  unsigned long ExtIter = config->GetExtIter();
  bool cont_adjoint     = config->GetContinuous_Adjoint();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_Flow() || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == ROE));
  bool limiter          = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool center           = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst       = center && (config->GetKind_Centered_Flow() == JST);
  bool engine           = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk    = ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0));
  bool nearfield        = (config->GetnMarker_NearFieldBound() != 0);
  bool interface        = (config->GetnMarker_InterfaceBound() != 0);
  bool fixed_cl         = config->GetFixed_CL_Mode();
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
  unsigned short kind_row_dissipation = config->GetKind_RoeLowDiss();
  bool roe_low_dissipation  = (kind_row_dissipation != NO_ROELOWDISS) && (config->GetKind_Upwind_Flow() == ROE);

  /*--- Update the angle of attack at the far-field for fixed CL calculations. ---*/
  
  if (fixed_cl) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }

  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

  /*--- Compute the engine properties ---*/

  if (engine) { GetPower_Properties(geometry, config, iMesh, Output); }

  /*--- Compute the actuator disk properties and distortion levels ---*/

  if (actuator_disk) {
    Set_MPI_ActDisk(solver_container, geometry, config);
    GetPower_Properties(geometry, config, iMesh, Output);
    SetActDisk_BCThrust(geometry, solver_container, config, iMesh, Output);
  }

  /*--- Compute Interface MPI ---*/

  if (interface) { Set_MPI_Interface(geometry, config); }

  /*--- Compute NearField MPI ---*/

  if (nearfield) { Set_MPI_Nearfield(geometry, config); }

 
  /*--- Upwind second order reconstruction ---*/
  
  if ((muscl && !center) && (iMesh == MESH_0) && !Output) {
    
    /*--- Gradient computation ---*/
    
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimitive_Gradient_GG(geometry, config);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimitive_Gradient_LS(geometry, config);
    }
    
    /*--- Limiter computation ---*/
    
    if (limiter && (iMesh == MESH_0)
        && !Output && !van_albada) { SetPrimitive_Limiter(geometry, config); }
    
  }
  
  /*--- Artificial dissipation ---*/
  
  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
  /*--- Roe Low Dissipation Sensor ---*/
  
  if (roe_low_dissipation){
    SetRoe_Dissipation(geometry, config);
    if (kind_row_dissipation == FD_DUCROS || kind_row_dissipation == NTS_DUCROS){
      SetUpwind_Ducros_Sensor(geometry, config);
    }
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (implicit && !disc_adjoint) Jacobian.SetValZero();

  /*--- Error message ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }
  
}

void CEulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh) { }

unsigned long CEulerSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the non-physical points vector ---*/
    
    node[iPoint]->SetNon_Physical(false);
    
    /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/
    
    RightSol = node[iPoint]->SetPrimVar(FluidModel);
    node[iPoint]->SetSecondaryVar(FluidModel);

    if (!RightSol) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  return ErrorCounter;
}
void CEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {
  
  su2double *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time,
  Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, ProjVel_i, ProjVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
    bool time_steping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetMax_Lambda_Inv(0.0);
  
  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    
    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
      
    }
  }
  
  /*--- Each element uses their own speed, steady state simulation ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
    
  }
  
  
  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;
    
    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
#endif
  }
  
  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  
  if (time_steping) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            /*--- Sets the regular CFL equal to the unsteady CFL ---*/
            config->SetCFL(iMesh,config->GetUnst_CFL());
            
            /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
             it computes the time step based on the unsteady CFL ---*/
            if (config->GetCFL(iMesh) == 0.0) {
                node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
            } else {
                node[iPoint]->SetDelta_Time(Global_Delta_Time);
            }
        }
  }
  
  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
    
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = rbuf_time;
#endif
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }
  
  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
  
  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  
}

void CEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iEdge, iPoint, jPoint;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool jst_scheme = ((config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0));
  bool grid_movement = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge, set normal vectors, and number of neighbors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
    
    /*--- Set primitive variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
    
    /*--- Set the largest convective eigenvalue ---*/
    
    numerics->SetLambda(node[iPoint]->GetLambda(), node[jPoint]->GetLambda());
    
    /*--- Set undivided laplacian an pressure based sensor ---*/
    
    if (jst_scheme) {
      numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
      numerics->SetSensor(node[iPoint]->GetSensor(), node[jPoint]->GetSensor());
    }
    
    /*--- Grid movement ---*/
    
    if (grid_movement) {
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    }
    
    /*--- Compute residuals, and Jacobians ---*/
    
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    
    /*--- Update convective and artificial dissipation residuals ---*/
    
    LinSysRes.AddBlock(iPoint, Res_Conv);
    LinSysRes.SubtractBlock(jPoint, Res_Conv);
    
    /*--- Set implicit computation ---*/
    if (implicit) {
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    }
  }
  
}

void CEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  su2double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, RoeVelocity[3] = {0.0,0.0,0.0}, R, sq_vel, RoeEnthalpy,
  *V_i, *V_j, *S_i, *S_j, *Limiter_i = NULL, *Limiter_j = NULL, sqvel, Non_Physical = 1.0, Sensor_i, Sensor_j, Dissipation_i, Dissipation_j, *Coord_i, *Coord_j;
  
  su2double z, velocity2_i, velocity2_j, mach_i, mach_j, vel_i_corr[3], vel_j_corr[3];
  
  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar;
  
  bool neg_density_i = false, neg_density_j = false, neg_pressure_i = false, neg_pressure_j = false, neg_sound_speed = false;
  
  unsigned long ExtIter = config->GetExtIter();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_Flow() && (iMesh == MESH_0));
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter          = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool grid_movement    = config->GetGrid_Movement();
  bool roe_turkel       = (config->GetKind_Upwind_Flow() == TURKEL);
  bool ideal_gas        = (config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS );
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
  bool low_mach_corr    = config->Low_Mach_Correction();
  unsigned short kind_dissipation = config->GetKind_RoeLowDiss();
    
  /*--- Loop over all the edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Roe Turkel preconditioning ---*/
    
    if (roe_turkel) {
      sqvel = 0.0;
      for (iDim = 0; iDim < nDim; iDim ++)
        sqvel += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
      numerics->SetVelocity2_Inf(sqvel);
    }
    
    /*--- Grid movement ---*/
    
    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    /*--- Get primitive variables ---*/
    
    V_i = node[iPoint]->GetPrimitive(); V_j = node[jPoint]->GetPrimitive();
    S_i = node[iPoint]->GetSecondary(); S_j = node[jPoint]->GetSecondary();

    /*--- High order reconstruction using MUSCL strategy ---*/
    
    if (muscl) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      if (limiter) {
        Limiter_i = node[iPoint]->GetLimiter_Primitive();
        Limiter_j = node[jPoint]->GetLimiter_Primitive();
      }
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        Non_Physical = node[iPoint]->GetNon_Physical()*node[jPoint]->GetNon_Physical();
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim]*Non_Physical;
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim]*Non_Physical;
        }
        if (limiter) {
          if (van_albada){
            Limiter_i[iVar] = (V_j[iVar]-V_i[iVar])*(2.0*Project_Grad_i + V_j[iVar]-V_i[iVar])/(4*Project_Grad_i*Project_Grad_i+(V_j[iVar]-V_i[iVar])*(V_j[iVar]-V_i[iVar])+EPS);
            Limiter_j[iVar] = (V_j[iVar]-V_i[iVar])*(-2.0*Project_Grad_j + V_j[iVar]-V_i[iVar])/(4*Project_Grad_j*Project_Grad_j+(V_j[iVar]-V_i[iVar])*(V_j[iVar]-V_i[iVar])+EPS);
          }
          Primitive_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Primitive_i[iVar] = V_i[iVar] + Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }

      /*--- Recompute the extrapolated quantities in a
       thermodynamic consistent way  ---*/

      if (!ideal_gas || low_mach_corr) { ComputeConsExtrapolation(config); }

      /*--- Low-Mach number correction ---*/

      if (low_mach_corr) {

        velocity2_i = 0.0;
        velocity2_j = 0.0;
        
        for (iDim = 0; iDim < nDim; iDim++) {
          velocity2_i += Primitive_i[iDim+1]*Primitive_i[iDim+1];
          velocity2_j += Primitive_j[iDim+1]*Primitive_j[iDim+1];
        }
        mach_i = sqrt(velocity2_i)/Primitive_i[nDim+4];
        mach_j = sqrt(velocity2_j)/Primitive_j[nDim+4];

        z = min(max(mach_i,mach_j),1.0);
        velocity2_i = 0.0;
        velocity2_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          vel_i_corr[iDim] = ( Primitive_i[iDim+1] + Primitive_j[iDim+1] )/2.0 \
                  + z * ( Primitive_i[iDim+1] - Primitive_j[iDim+1] )/2.0;
          vel_j_corr[iDim] = ( Primitive_i[iDim+1] + Primitive_j[iDim+1] )/2.0 \
                  + z * ( Primitive_j[iDim+1] - Primitive_i[iDim+1] )/2.0;

          velocity2_j += vel_j_corr[iDim]*vel_j_corr[iDim];
          velocity2_i += vel_i_corr[iDim]*vel_i_corr[iDim];

          Primitive_i[iDim+1] = vel_i_corr[iDim];
          Primitive_j[iDim+1] = vel_j_corr[iDim];
        }

        FluidModel->SetEnergy_Prho(Primitive_i[nDim+1],Primitive_i[nDim+2]);
        Primitive_i[nDim+3]= FluidModel->GetStaticEnergy() + Primitive_i[nDim+1]/Primitive_i[nDim+2] + 0.5*velocity2_i;
        
        FluidModel->SetEnergy_Prho(Primitive_j[nDim+1],Primitive_j[nDim+2]);
        Primitive_j[nDim+3]= FluidModel->GetStaticEnergy() + Primitive_j[nDim+1]/Primitive_j[nDim+2] + 0.5*velocity2_j;
        
      }
      
      /*--- Check for non-physical solutions after reconstruction. If found,
       use the cell-average value of the solution. This results in a locally
       first-order approximation, but this is typically only active
       during the start-up of a calculation. If non-physical, use the 
       cell-averaged state. ---*/
      
      neg_pressure_i = (Primitive_i[nDim+1] < 0.0); neg_pressure_j = (Primitive_j[nDim+1] < 0.0);
      neg_density_i  = (Primitive_i[nDim+2] < 0.0); neg_density_j  = (Primitive_j[nDim+2] < 0.0);

      R = sqrt(fabs(Primitive_j[nDim+2]/Primitive_i[nDim+2]));
      sq_vel = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        RoeVelocity[iDim] = (R*Primitive_j[iDim+1]+Primitive_i[iDim+1])/(R+1);
        sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
      }
      RoeEnthalpy = (R*Primitive_j[nDim+3]+Primitive_i[nDim+3])/(R+1);
      neg_sound_speed = ((Gamma-1)*(RoeEnthalpy-0.5*sq_vel) < 0.0);
      
      if (neg_sound_speed) {
        for (iVar = 0; iVar < nPrimVar; iVar++) {
          Primitive_i[iVar] = V_i[iVar];
          Primitive_j[iVar] = V_j[iVar]; }
        Secondary_i[0] = S_i[0]; Secondary_i[1] = S_i[1];
        Secondary_j[0] = S_i[0]; Secondary_j[1] = S_i[1];
        counter_local++;
      }
      
      if (neg_density_i || neg_pressure_i) {
        for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = V_i[iVar];
        Secondary_i[0] = S_i[0]; Secondary_i[1] = S_i[1];
        counter_local++;
      }
      
      if (neg_density_j || neg_pressure_j) {
        for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = V_j[iVar];
        Secondary_j[0] = S_j[0]; Secondary_j[1] = S_j[1];
        counter_local++;
      }

      numerics->SetPrimitive(Primitive_i, Primitive_j);
      numerics->SetSecondary(Secondary_i, Secondary_j);
      
    }
    else {
      
      /*--- Set conservative variables without reconstruction ---*/
      
      numerics->SetPrimitive(V_i, V_j);
      numerics->SetSecondary(S_i, S_j);
      
    }
    
    /*--- Roe Low Dissipation Scheme ---*/
    
    if (kind_dissipation != NO_ROELOWDISS){
      
      Dissipation_i = node[iPoint]->GetRoe_Dissipation();
      Dissipation_j = node[jPoint]->GetRoe_Dissipation();
      numerics->SetDissipation(Dissipation_i, Dissipation_j);
            
      if (kind_dissipation == FD_DUCROS || kind_dissipation == NTS_DUCROS){
        Sensor_i = node[iPoint]->GetSensor();
        Sensor_j = node[jPoint]->GetSensor();
        numerics->SetSensor(Sensor_i, Sensor_j);
      }
      if (kind_dissipation == NTS || kind_dissipation == NTS_DUCROS){
        Coord_i = geometry->node[iPoint]->GetCoord();
        Coord_j = geometry->node[jPoint]->GetCoord();
        numerics->SetCoord(Coord_i, Coord_j);
      }
    }
      
    /*--- Compute the residual ---*/
    
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

    /*--- Update residual value ---*/
    
    LinSysRes.AddBlock(iPoint, Res_Conv);
    LinSysRes.SubtractBlock(jPoint, Res_Conv);
    
    /*--- Set implicit Jacobians ---*/
    
    if (implicit) {
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    }
    
    /*--- Roe Turkel preconditioning, set the value of beta ---*/
    
    if (roe_turkel) {
      node[iPoint]->SetPreconditioner_Beta(numerics->GetPrecond_Beta());
      node[jPoint]->SetPreconditioner_Beta(numerics->GetPrecond_Beta());
    }
    
    /*--- Set the final value of the Roe dissipation coefficient ---*/
    
    if (kind_dissipation != NO_ROELOWDISS){
      node[iPoint]->SetRoe_Dissipation(numerics->GetDissipation());
      node[jPoint]->SetRoe_Dissipation(numerics->GetDissipation());      
    }
    
  }

  /*--- Warning message about non-physical reconstructions ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Reconstr(counter_global);
  }
}

void CEulerSolver::ComputeConsExtrapolation(CConfig *config) {
  
  unsigned short iDim;
  
  su2double density_i = Primitive_i[nDim+2];
  su2double pressure_i = Primitive_i[nDim+1];
  su2double velocity2_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    velocity2_i += Primitive_i[iDim+1]*Primitive_i[iDim+1];
  }
  
  FluidModel->SetTDState_Prho(pressure_i, density_i);
  
  Primitive_i[0]= FluidModel->GetTemperature();
  Primitive_i[nDim+3]= FluidModel->GetStaticEnergy() + Primitive_i[nDim+1]/Primitive_i[nDim+2] + 0.5*velocity2_i;
  Primitive_i[nDim+4]= FluidModel->GetSoundSpeed();
  Secondary_i[0]=FluidModel->GetdPdrho_e();
  Secondary_i[1]=FluidModel->GetdPde_rho();
  
  
  su2double density_j = Primitive_j[nDim+2];
  su2double pressure_j = Primitive_j[nDim+1];
  su2double velocity2_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    velocity2_j += Primitive_j[iDim+1]*Primitive_j[iDim+1];
  }
  
  FluidModel->SetTDState_Prho(pressure_j, density_j);
  
  Primitive_j[0]= FluidModel->GetTemperature();
  Primitive_j[nDim+3]= FluidModel->GetStaticEnergy() + Primitive_j[nDim+1]/Primitive_j[nDim+2] + 0.5*velocity2_j;
  Primitive_j[nDim+4]=FluidModel->GetSoundSpeed();
  Secondary_j[0]=FluidModel->GetdPdrho_e();
  Secondary_j[1]=FluidModel->GetdPde_rho();
  
}

void CEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint;
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool rotating_frame   = config->GetRotating_Frame();
  bool axisymmetric     = config->GetAxisymmetric();
  bool gravity          = (config->GetGravityForce() == YES);
  bool harmonic_balance = (config->GetUnsteady_Simulation() == HARMONIC_BALANCE);
  bool windgust         = config->GetWind_Gust();
  bool body_force       = config->GetBody_Force();

  /*--- Initialize the source residual to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;

  if (body_force) {

    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the conservative variables ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(),
                                node[iPoint]->GetSolution());

      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, config);

      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
  }

  if (rotating_frame) {

    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the conservative variables ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(),
                                node[iPoint]->GetSolution());
      
      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  if (axisymmetric) {
    
    /*--- Zero out Jacobian structure ---*/
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Set solution  ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Set y coordinate ---*/
      numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
      
      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  }
  
  if (gravity) {
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Set solution  ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, config);
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
    
  }
  
  if (harmonic_balance) {
    
    su2double Volume, Source;
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Get control volume ---*/
      Volume = geometry->node[iPoint]->GetVolume();
      
      /*--- Get stored time spectral source term ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Source = node[iPoint]->GetHarmonicBalance_Source(iVar);        
        Residual[iVar] = Source*Volume;
      }
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
  }
  
  if (windgust) {
    
    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the wind gust ---*/
      numerics->SetWindGust(node[iPoint]->GetWindGust(), node[iPoint]->GetWindGust());
      
      /*--- Load the wind gust derivatives ---*/
      numerics->SetWindGustDer(node[iPoint]->GetWindGustDer(), node[iPoint]->GetWindGustDer());
      
      /*--- Load the primitive variables ---*/
      numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[iPoint]->GetPrimitive());
      
      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
}

void CEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  /* This method should be used to call any new source terms for a particular problem*/
  /* This method calls the new child class in CNumerics, where the new source term should be implemented.  */
  
  /* Next we describe how to get access to some important quanties for this method */
  /* Access to all points in the current geometric mesh by saying: nPointDomain */
  /* Get the vector of conservative variables at some point iPoint = node[iPoint]->GetSolution() */
  /* Get the volume (or area in 2D) associated with iPoint = node[iPoint]->GetVolume() */
  /* Get the vector of geometric coordinates of point iPoint = node[iPoint]->GetCoord() */
  
}

void CEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {
  
  su2double *Normal, Area, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda,
  ProjVel, ProjVel_i, ProjVel_j, *GridVel, *GridVel_i, *GridVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;
  
  bool grid_movement = config->GetGrid_Movement();
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetLambda(0.0);
  }
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;

    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddLambda(Lambda);
      }
      
    }
  }
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_MaxEigenvalue(geometry, config);
  
}

void CEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint, jPoint, iEdge;
  su2double Pressure_i = 0, Pressure_j = 0, *Diff;
  unsigned short iVar;
  bool boundary_i, boundary_j;
    
  Diff = new su2double[nVar];
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetUnd_LaplZero();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Solution differences ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
    
    /*--- Correction for compressible flows which use the enthalpy ---*/
    
    Pressure_i = node[iPoint]->GetPressure();
    Pressure_j = node[jPoint]->GetPressure();
    Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1) + Pressure_i) - (node[jPoint]->GetSolution(nVar-1) + Pressure_j);
    
    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both in the boundary ---*/
    
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    
  }
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_Undivided_Laplacian(geometry, config);
  
  delete [] Diff;
  
}

void CEulerSolver::SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  su2double Pressure_i = 0.0, Pressure_j = 0.0;
  bool boundary_i, boundary_j;
  
  /*--- Reset variables to store the undivided pressure ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    iPoint_UndLapl[iPoint] = 0.0;
    jPoint_UndLapl[iPoint] = 0.0;
  }
  
  /*--- Evaluate the pressure sensor ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Pressure_i = node[iPoint]->GetPressure();
    Pressure_j = node[jPoint]->GetPressure();
    
    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both on the boundary ---*/
    
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }
    
  }
  
  /*--- Set pressure switch for each point ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetSensor(fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_Sensor(geometry, config);
  
}

void CEulerSolver::SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config){
  
  unsigned long iPoint, jPoint;
  unsigned short iNeigh, iDim;
  
  su2double *Vorticity;
  
  su2double uixi = 0.0, Ducros_i = 0.0, Ducros_j = 0.0, Omega = 0.0;
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
    
    /*---- Dilatation for iPoint ---*/
    
    uixi=0.0;
    for(iDim = 0; iDim < nDim; iDim++){
      uixi += node[iPoint]->GetGradient_Primitive(iDim+1, iDim);
    } 
    
    /*--- Compute norm of vorticity ---*/
    
    Vorticity = node[iPoint]->GetVorticity();    
    Omega = 0.0;
    for (iDim = 0; iDim < nDim; iDim++){
      Omega += Vorticity[iDim]*Vorticity[iDim];
    }
    Omega = sqrt(Omega);
    
    /*---- Ducros sensor for iPoint ---*/
    
    if (config->GetKind_RoeLowDiss() == FD_DUCROS){
      Ducros_i = -uixi / (fabs(uixi) + Omega + 1e-20);
    } else if (config->GetKind_RoeLowDiss() == NTS_DUCROS){
      Ducros_i = pow(uixi,2.0) /(pow(uixi,2.0)+ pow(Omega,2.0) + 1e-20);
    }
    
    node[iPoint]->SetSensor(Ducros_i);
    
    /*---- Ducros sensor for neighbor points of iPoint to avoid lower the dissipation in regions near the shock ---*/
    
    for (iNeigh = 0; iNeigh > geometry->node[iPoint]->GetnNeighbor(); iNeigh++){
      
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      
      /*---- Dilatation for jPoint ---*/
      
      uixi=0.0;
      for(iDim = 0; iDim < nDim; iDim++){
        uixi += node[jPoint]->GetGradient_Primitive(iDim+1, iDim);
      } 
      
      /*--- Compute norm of vorticity ---*/
      
      Vorticity = node[jPoint]->GetVorticity();      
      Omega = 0.0;
      for (iDim = 0; iDim < nDim; iDim++){
        Omega += Vorticity[iDim]*Vorticity[iDim];
      }
      Omega = sqrt(Omega);
      
      if (config->GetKind_RoeLowDiss() == FD_DUCROS){
        Ducros_j = -uixi / (fabs(uixi) + Omega + 1e-20);
      } else if (config->GetKind_RoeLowDiss() == NTS_DUCROS){
        Ducros_j = pow(uixi,2.0) /(pow(uixi,2.0)+ pow(Omega,2.0) + 1e-20);
      }      
      node[iPoint]->SetSensor(max(node[iPoint]->GetSensor(), Ducros_j));

    }
  }
  
  Set_MPI_Sensor(geometry, config);
  
}

void CEulerSolver::Pressure_Forces(CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double Pressure = 0.0, *Normal = NULL, MomentDist[3] = {0.0,0.0,0.0}, *Coord,
  factor, NFPressOF, RefVel2, RefTemp, RefDensity, RefPressure, Mach2Vel, Mach_Motion,
  Force[3] = {0.0,0.0,0.0};
  string Marker_Tag, Monitoring_Tag;
  su2double MomentX_Force[3] = {0.0,0.0,0.0}, MomentY_Force[3] = {0.0,0.0,0.0}, MomentZ_Force[3] = {0.0,0.0,0.0};
  su2double AxiFactor;

#ifdef HAVE_MPI
  su2double MyAllBound_CD_Inv, MyAllBound_CL_Inv, MyAllBound_CSF_Inv, MyAllBound_CMx_Inv, MyAllBound_CMy_Inv, MyAllBound_CMz_Inv, MyAllBound_CoPx_Inv, MyAllBound_CoPy_Inv, MyAllBound_CoPz_Inv, MyAllBound_CFx_Inv, MyAllBound_CFy_Inv, MyAllBound_CFz_Inv, MyAllBound_CT_Inv, MyAllBound_CQ_Inv, MyAllBound_CNearFieldOF_Inv, *MySurface_CL_Inv = NULL, *MySurface_CD_Inv = NULL, *MySurface_CSF_Inv = NULL, *MySurface_CEff_Inv = NULL, *MySurface_CFx_Inv = NULL, *MySurface_CFy_Inv = NULL, *MySurface_CFz_Inv = NULL, *MySurface_CMx_Inv = NULL, *MySurface_CMy_Inv = NULL, *MySurface_CMz_Inv = NULL;
#endif
  
  su2double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta            = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea    = config->GetRefArea();
  su2double RefLength = config->GetRefLength();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  su2double *Origin = NULL;
  if (config->GetnMarker_Monitoring() != 0){
    Origin = config->GetRefOriginMoment(0);
  }
  bool grid_movement        = config->GetGrid_Movement();
  bool axisymmetric         = config->GetAxisymmetric();

  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream
   values, which is the standard convention. ---*/
  
  RefTemp     = Temperature_Inf;
  RefDensity  = Density_Inf;
  RefPressure = Pressure_Inf;
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);
  
  /*-- Variables initialization ---*/
  
  Total_CD = 0.0;           Total_CL = 0.0;    Total_CSF = 0.0;     Total_CEff = 0.0;
  Total_CMx = 0.0;          Total_CMy = 0.0;   Total_CMz = 0.0;
  Total_CoPx = 0.0;         Total_CoPy = 0.0;  Total_CoPz = 0.0;
  Total_CFx = 0.0;          Total_CFy = 0.0;   Total_CFz = 0.0;
  Total_CT = 0.0;           Total_CQ = 0.0;    Total_CMerit = 0.0;
  Total_CNearFieldOF = 0.0; Total_Heat = 0.0;  Total_MaxHeat = 0.0;
  
  AllBound_CD_Inv = 0.0;        AllBound_CL_Inv = 0.0; AllBound_CSF_Inv = 0.0;
  AllBound_CMx_Inv = 0.0;          AllBound_CMy_Inv = 0.0;   AllBound_CMz_Inv = 0.0;
  AllBound_CoPx_Inv = 0.0;         AllBound_CoPy_Inv = 0.0;  AllBound_CoPz_Inv = 0.0;
  AllBound_CFx_Inv = 0.0;          AllBound_CFy_Inv = 0.0;   AllBound_CFz_Inv = 0.0;
  AllBound_CT_Inv = 0.0;           AllBound_CQ_Inv = 0.0;    AllBound_CMerit_Inv = 0.0;
  AllBound_CNearFieldOF_Inv = 0.0; AllBound_CEff_Inv = 0.0;
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Inv[iMarker_Monitoring]      = 0.0; Surface_CD_Inv[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring] = 0.0; Surface_CEff_Inv[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring]        = 0.0; Surface_CFy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring]        = 0.0; Surface_CMx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]        = 0.0; Surface_CMz_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CL[iMarker_Monitoring]          = 0.0; Surface_CD[iMarker_Monitoring]          = 0.0;
    Surface_CSF[iMarker_Monitoring]     = 0.0; Surface_CEff[iMarker_Monitoring]           = 0.0;
    Surface_CFx[iMarker_Monitoring]            = 0.0; Surface_CFy[iMarker_Monitoring]            = 0.0;
    Surface_CFz[iMarker_Monitoring]            = 0.0; Surface_CMx[iMarker_Monitoring]            = 0.0;
    Surface_CMy[iMarker_Monitoring]            = 0.0; Surface_CMz[iMarker_Monitoring]            = 0.0;
  }
  
  /*--- Loop over the Euler and Navier-Stokes markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    /*--- Obtain the origin for the moment computation for a particular marker ---*/
    
    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }
    
    if ((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
        (Boundary == ISOTHERMAL) || (Boundary == NEARFIELD_BOUNDARY) ||
        (Boundary == INLET_FLOW) || (Boundary == OUTLET_FLOW) ||
        (Boundary == ACTDISK_INLET) || (Boundary == ACTDISK_OUTLET)||
        (Boundary == ENGINE_INFLOW) || (Boundary == ENGINE_EXHAUST)) {
      
      /*--- Forces initialization at each Marker ---*/
      
      CD_Inv[iMarker] = 0.0;        CL_Inv[iMarker] = 0.0; CSF_Inv[iMarker] = 0.0;
      CMx_Inv[iMarker] = 0.0;          CMy_Inv[iMarker] = 0.0;   CMz_Inv[iMarker] = 0.0;
      CoPx_Inv[iMarker] = 0.0;          CoPy_Inv[iMarker] = 0.0;   CoPz_Inv[iMarker] = 0.0;
      CFx_Inv[iMarker] = 0.0;          CFy_Inv[iMarker] = 0.0;   CFz_Inv[iMarker] = 0.0;
      CT_Inv[iMarker] = 0.0;           CQ_Inv[iMarker] = 0.0;    CMerit_Inv[iMarker] = 0.0;
      CNearFieldOF_Inv[iMarker] = 0.0; CEff_Inv[iMarker] = 0.0;
      
      for (iDim = 0; iDim < nDim; iDim++) ForceInviscid[iDim] = 0.0;
      MomentInviscid[0] = 0.0; MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;
      MomentX_Force[0] = 0.0; MomentX_Force[1] = 0.0; MomentX_Force[2] = 0.0;
      MomentY_Force[0] = 0.0; MomentY_Force[1] = 0.0; MomentY_Force[2] = 0.0;
      MomentZ_Force[0] = 0.0; MomentZ_Force[1] = 0.0; MomentZ_Force[2] = 0.0;

      NFPressOF = 0.0;
      
      /*--- Loop over the vertices to compute the forces ---*/
      
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        Pressure = node[iPoint]->GetPressure();
        
        CPressure[iMarker][iVertex] = (Pressure - RefPressure)*factor*RefArea;
        
        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/
        
        if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {
          
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          
          /*--- Quadratic objective function for the near-field.
           This uses the infinity pressure regardless of Mach number. ---*/
          
          NFPressOF += 0.5*(Pressure - Pressure_Inf)*(Pressure - Pressure_Inf)*Normal[nDim-1];
          
          for (iDim = 0; iDim < nDim; iDim++) {
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }
          
          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = -(Pressure - Pressure_Inf) * Normal[iDim] * factor * AxiFactor;
            ForceInviscid[iDim] += Force[iDim];
          }
          
          /*--- Moment with respect to the reference axis ---*/
          
          if (nDim == 3) {
            MomentInviscid[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLength;
            MomentX_Force[1]  += (-Force[1]*Coord[2]);
            MomentX_Force[2]  += (Force[2]*Coord[1]);

            MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLength;
            MomentY_Force[2]  += (-Force[2]*Coord[0]);
            MomentY_Force[0]  += (Force[0]*Coord[2]);
          }
          MomentInviscid[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLength;
          MomentZ_Force[0]  += (-Force[0]*Coord[1]);
          MomentZ_Force[1]  += (Force[1]*Coord[0]);
        }
        
      }
      
      /*--- Project forces and store the non-dimensional coefficients ---*/
      
      if (Monitoring == YES) {
        
        if (Boundary != NEARFIELD_BOUNDARY) {
          if (nDim == 2) {
            CD_Inv[iMarker]     =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
            CL_Inv[iMarker]     = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
            CEff_Inv[iMarker]   = CL_Inv[iMarker] / (CD_Inv[iMarker]+EPS);
            CMz_Inv[iMarker]    = MomentInviscid[2];
            CoPx_Inv[iMarker]   = MomentZ_Force[1];
            CoPy_Inv[iMarker]   = -MomentZ_Force[0];
            CFx_Inv[iMarker]    = ForceInviscid[0];
            CFy_Inv[iMarker]    = ForceInviscid[1];
            CT_Inv[iMarker]     = -CFx_Inv[iMarker];
            CQ_Inv[iMarker]     = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker] = CT_Inv[iMarker] / (CQ_Inv[iMarker] + EPS);
          }
          if (nDim == 3) {
            CD_Inv[iMarker]      =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
            CL_Inv[iMarker]      = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
            CSF_Inv[iMarker]     = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
            CEff_Inv[iMarker]    = CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
            CMx_Inv[iMarker]     = MomentInviscid[0];
            CMy_Inv[iMarker]     = MomentInviscid[1];
            CMz_Inv[iMarker]     = MomentInviscid[2];
            CoPx_Inv[iMarker]    = -MomentY_Force[0];
            CoPz_Inv[iMarker]    = MomentY_Force[2];
            CFx_Inv[iMarker]     = ForceInviscid[0];
            CFy_Inv[iMarker]     = ForceInviscid[1];
            CFz_Inv[iMarker]     = ForceInviscid[2];
            CT_Inv[iMarker]      = -CFz_Inv[iMarker];
            CQ_Inv[iMarker]      = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker]  = CT_Inv[iMarker] / (CQ_Inv[iMarker] + EPS);
          }
          
          AllBound_CD_Inv           += CD_Inv[iMarker];
          AllBound_CL_Inv           += CL_Inv[iMarker];
          AllBound_CSF_Inv          += CSF_Inv[iMarker];
          AllBound_CEff_Inv          = AllBound_CL_Inv / (AllBound_CD_Inv + EPS);
          AllBound_CMx_Inv          += CMx_Inv[iMarker];
          AllBound_CMy_Inv          += CMy_Inv[iMarker];
          AllBound_CMz_Inv          += CMz_Inv[iMarker];
          AllBound_CoPx_Inv         += CoPx_Inv[iMarker];
          AllBound_CoPy_Inv         += CoPy_Inv[iMarker];
          AllBound_CoPz_Inv         += CoPz_Inv[iMarker];
          AllBound_CFx_Inv          += CFx_Inv[iMarker];
          AllBound_CFy_Inv          += CFy_Inv[iMarker];
          AllBound_CFz_Inv          += CFz_Inv[iMarker];
          AllBound_CT_Inv           += CT_Inv[iMarker];
          AllBound_CQ_Inv           += CQ_Inv[iMarker];
          AllBound_CMerit_Inv        = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);
          
          /*--- Compute the coefficients per surface ---*/
          
          for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
            Marker_Tag = config->GetMarker_All_TagBound(iMarker);
            if (Marker_Tag == Monitoring_Tag) {
              Surface_CL_Inv[iMarker_Monitoring]      += CL_Inv[iMarker];
              Surface_CD_Inv[iMarker_Monitoring]      += CD_Inv[iMarker];
              Surface_CSF_Inv[iMarker_Monitoring]     += CSF_Inv[iMarker];
              Surface_CEff_Inv[iMarker_Monitoring]     = CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
              Surface_CFx_Inv[iMarker_Monitoring]     += CFx_Inv[iMarker];
              Surface_CFy_Inv[iMarker_Monitoring]     += CFy_Inv[iMarker];
              Surface_CFz_Inv[iMarker_Monitoring]     += CFz_Inv[iMarker];
              Surface_CMx_Inv[iMarker_Monitoring]     += CMx_Inv[iMarker];
              Surface_CMy_Inv[iMarker_Monitoring]     += CMy_Inv[iMarker];
              Surface_CMz_Inv[iMarker_Monitoring]     += CMz_Inv[iMarker];
            }
          }
          
        }
        
        /*--- At the Nearfield SU2 only cares about the pressure coeffient ---*/
        
        else {
          CNearFieldOF_Inv[iMarker] = NFPressOF;
          AllBound_CNearFieldOF_Inv += CNearFieldOF_Inv[iMarker];
        }
        
      }
      
      
    }
  }
  
#ifdef HAVE_MPI
  
  /*--- Add AllBound information using all the nodes ---*/
  
  MyAllBound_CD_Inv        = AllBound_CD_Inv;        AllBound_CD_Inv = 0.0;
  MyAllBound_CL_Inv        = AllBound_CL_Inv;        AllBound_CL_Inv = 0.0;
  MyAllBound_CSF_Inv   = AllBound_CSF_Inv;   AllBound_CSF_Inv = 0.0;
  AllBound_CEff_Inv = 0.0;
  MyAllBound_CMx_Inv          = AllBound_CMx_Inv;          AllBound_CMx_Inv = 0.0;
  MyAllBound_CMy_Inv          = AllBound_CMy_Inv;          AllBound_CMy_Inv = 0.0;
  MyAllBound_CMz_Inv          = AllBound_CMz_Inv;          AllBound_CMz_Inv = 0.0;
  MyAllBound_CoPx_Inv          = AllBound_CoPx_Inv;          AllBound_CoPx_Inv = 0.0;
  MyAllBound_CoPy_Inv          = AllBound_CoPy_Inv;          AllBound_CoPy_Inv = 0.0;
  MyAllBound_CoPz_Inv          = AllBound_CoPz_Inv;          AllBound_CoPz_Inv = 0.0;
  MyAllBound_CFx_Inv          = AllBound_CFx_Inv;          AllBound_CFx_Inv = 0.0;
  MyAllBound_CFy_Inv          = AllBound_CFy_Inv;          AllBound_CFy_Inv = 0.0;
  MyAllBound_CFz_Inv          = AllBound_CFz_Inv;          AllBound_CFz_Inv = 0.0;
  MyAllBound_CT_Inv           = AllBound_CT_Inv;           AllBound_CT_Inv = 0.0;
  MyAllBound_CQ_Inv           = AllBound_CQ_Inv;           AllBound_CQ_Inv = 0.0;
  AllBound_CMerit_Inv = 0.0;
  MyAllBound_CNearFieldOF_Inv = AllBound_CNearFieldOF_Inv; AllBound_CNearFieldOF_Inv = 0.0;
  
  SU2_MPI::Allreduce(&MyAllBound_CD_Inv, &AllBound_CD_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CL_Inv, &AllBound_CL_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CSF_Inv, &AllBound_CSF_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Inv = AllBound_CL_Inv / (AllBound_CD_Inv + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CMx_Inv, &AllBound_CMx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMy_Inv, &AllBound_CMy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMz_Inv, &AllBound_CMz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CoPx_Inv, &AllBound_CoPx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CoPy_Inv, &AllBound_CoPy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CoPz_Inv, &AllBound_CoPz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFx_Inv, &AllBound_CFx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFy_Inv, &AllBound_CFy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFz_Inv, &AllBound_CFz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CT_Inv, &AllBound_CT_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CQ_Inv, &AllBound_CQ_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CMerit_Inv = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CNearFieldOF_Inv, &AllBound_CNearFieldOF_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  /*--- Add the forces on the surfaces using all the nodes ---*/
  
  MySurface_CL_Inv      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Inv      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Inv = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_CL_Inv[iMarker_Monitoring]      = Surface_CL_Inv[iMarker_Monitoring];
    MySurface_CD_Inv[iMarker_Monitoring]      = Surface_CD_Inv[iMarker_Monitoring];
    MySurface_CSF_Inv[iMarker_Monitoring] = Surface_CSF_Inv[iMarker_Monitoring];
    MySurface_CEff_Inv[iMarker_Monitoring]       = Surface_CEff_Inv[iMarker_Monitoring];
    MySurface_CFx_Inv[iMarker_Monitoring]        = Surface_CFx_Inv[iMarker_Monitoring];
    MySurface_CFy_Inv[iMarker_Monitoring]        = Surface_CFy_Inv[iMarker_Monitoring];
    MySurface_CFz_Inv[iMarker_Monitoring]        = Surface_CFz_Inv[iMarker_Monitoring];
    MySurface_CMx_Inv[iMarker_Monitoring]        = Surface_CMx_Inv[iMarker_Monitoring];
    MySurface_CMy_Inv[iMarker_Monitoring]        = Surface_CMy_Inv[iMarker_Monitoring];
    MySurface_CMz_Inv[iMarker_Monitoring]        = Surface_CMz_Inv[iMarker_Monitoring];

    Surface_CL_Inv[iMarker_Monitoring]         = 0.0;
    Surface_CD_Inv[iMarker_Monitoring]         = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CEff_Inv[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Inv[iMarker_Monitoring]        = 0.0;
  }
  
  SU2_MPI::Allreduce(MySurface_CL_Inv, Surface_CL_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CD_Inv, Surface_CD_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CSF_Inv, Surface_CSF_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
    Surface_CEff_Inv[iMarker_Monitoring] = Surface_CL_Inv[iMarker_Monitoring] / (Surface_CD_Inv[iMarker_Monitoring] + EPS);
  SU2_MPI::Allreduce(MySurface_CFx_Inv, Surface_CFx_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFy_Inv, Surface_CFy_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFz_Inv, Surface_CFz_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMx_Inv, Surface_CMx_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMy_Inv, Surface_CMy_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMz_Inv, Surface_CMz_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  delete [] MySurface_CL_Inv; delete [] MySurface_CD_Inv; delete [] MySurface_CSF_Inv;
  delete [] MySurface_CEff_Inv;  delete [] MySurface_CFx_Inv;   delete [] MySurface_CFy_Inv;
  delete [] MySurface_CFz_Inv;   delete [] MySurface_CMx_Inv;   delete [] MySurface_CMy_Inv;
  delete [] MySurface_CMz_Inv;
  
#endif
  
  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/
  
  Total_CD            = AllBound_CD_Inv;
  Total_CL            = AllBound_CL_Inv;
  Total_CSF           = AllBound_CSF_Inv;
  Total_CEff          = Total_CL / (Total_CD + EPS);
  Total_CFx           = AllBound_CFx_Inv;
  Total_CFy           = AllBound_CFy_Inv;
  Total_CFz           = AllBound_CFz_Inv;
  Total_CMx           = AllBound_CMx_Inv;
  Total_CMy           = AllBound_CMy_Inv;
  Total_CMz           = AllBound_CMz_Inv;
  Total_CoPx          = AllBound_CoPx_Inv;
  Total_CoPy          = AllBound_CoPy_Inv;
  Total_CoPz          = AllBound_CoPz_Inv;
  Total_CT            = AllBound_CT_Inv;
  Total_CQ            = AllBound_CQ_Inv;
  Total_CMerit        = Total_CT / (Total_CQ + EPS);
  Total_CNearFieldOF  = AllBound_CNearFieldOF_Inv;
  
  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]      = Surface_CL_Inv[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]      = Surface_CD_Inv[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring] = Surface_CSF_Inv[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       = Surface_CL_Inv[iMarker_Monitoring] / (Surface_CD_Inv[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        = Surface_CFx_Inv[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        = Surface_CFy_Inv[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        = Surface_CFz_Inv[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        = Surface_CMx_Inv[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        = Surface_CMy_Inv[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        = Surface_CMz_Inv[iMarker_Monitoring];
  }
  
}

void CEulerSolver::Momentum_Forces(CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double *Normal = NULL, MomentDist[3] = {0.0,0.0,0.0}, *Coord, Area,
  factor, RefVel2, RefTemp, RefDensity,  Mach2Vel, Mach_Motion,
  Force[3] = {0.0,0.0,0.0}, Velocity[3], MassFlow, Density;
  string Marker_Tag, Monitoring_Tag;
  su2double MomentX_Force[3] = {0.0,0.0,0.0}, MomentY_Force[3] = {0.0,0.0,0.0}, MomentZ_Force[3] = {0.0,0.0,0.0};
  su2double AxiFactor;
  
#ifdef HAVE_MPI
  su2double MyAllBound_CD_Mnt, MyAllBound_CL_Mnt, MyAllBound_CSF_Mnt,
  MyAllBound_CMx_Mnt, MyAllBound_CMy_Mnt, MyAllBound_CMz_Mnt,
  MyAllBound_CoPx_Mnt, MyAllBound_CoPy_Mnt, MyAllBound_CoPz_Mnt,
  MyAllBound_CFx_Mnt, MyAllBound_CFy_Mnt, MyAllBound_CFz_Mnt, MyAllBound_CT_Mnt,
  MyAllBound_CQ_Mnt,
  *MySurface_CL_Mnt = NULL, *MySurface_CD_Mnt = NULL, *MySurface_CSF_Mnt = NULL,
  *MySurface_CEff_Mnt = NULL, *MySurface_CFx_Mnt = NULL, *MySurface_CFy_Mnt = NULL,
  *MySurface_CFz_Mnt = NULL,
  *MySurface_CMx_Mnt = NULL, *MySurface_CMy_Mnt = NULL,  *MySurface_CMz_Mnt = NULL;
#endif
  
  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea     = config->GetRefArea();
  su2double RefLength  = config->GetRefLength();
  su2double Gas_Constant     = config->GetGas_ConstantND();
  su2double *Origin = NULL;
  if (config->GetnMarker_Monitoring() != 0){
    Origin = config->GetRefOriginMoment(0);
  }
  bool grid_movement         = config->GetGrid_Movement();
  bool axisymmetric          = config->GetAxisymmetric();
  
  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/
  
  RefTemp     = Temperature_Inf;
  RefDensity  = Density_Inf;
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);
  
  /*-- Variables initialization ---*/
  
  AllBound_CD_Mnt = 0.0;        AllBound_CL_Mnt = 0.0; AllBound_CSF_Mnt = 0.0;
  AllBound_CMx_Mnt = 0.0;          AllBound_CMy_Mnt = 0.0;   AllBound_CMz_Mnt = 0.0;
  AllBound_CoPx_Mnt = 0.0;          AllBound_CoPy_Mnt = 0.0;   AllBound_CoPz_Mnt = 0.0;
  AllBound_CFx_Mnt = 0.0;          AllBound_CFy_Mnt = 0.0;   AllBound_CFz_Mnt = 0.0;
  AllBound_CT_Mnt = 0.0;           AllBound_CQ_Mnt = 0.0;    AllBound_CMerit_Mnt = 0.0;
  AllBound_CEff_Mnt = 0.0;
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Mnt[iMarker_Monitoring]      = 0.0; Surface_CD_Mnt[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Mnt[iMarker_Monitoring] = 0.0; Surface_CEff_Mnt[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Mnt[iMarker_Monitoring]        = 0.0; Surface_CFy_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Mnt[iMarker_Monitoring]        = 0.0; Surface_CMy_Mnt[iMarker_Monitoring]        = 0.0; Surface_CMz_Mnt[iMarker_Monitoring]        = 0.0;
  }
  
  /*--- Loop over the Inlet -Outlet Markers  ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    /*--- Obtain the origin for the moment computation for a particular marker ---*/
    
    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }
    
    if ((Boundary == INLET_FLOW) || (Boundary == OUTLET_FLOW) ||
        (Boundary == ACTDISK_INLET) || (Boundary == ACTDISK_OUTLET)||
        (Boundary == ENGINE_INFLOW) || (Boundary == ENGINE_EXHAUST)) {
      
      /*--- Forces initialization at each Marker ---*/
      
      CD_Mnt[iMarker] = 0.0;           CL_Mnt[iMarker] = 0.0;    CSF_Mnt[iMarker] = 0.0;
      CFx_Mnt[iMarker] = 0.0;          CFy_Mnt[iMarker] = 0.0;   CFz_Mnt[iMarker] = 0.0;
      CMx_Mnt[iMarker] = 0.0;          CMy_Mnt[iMarker] = 0.0;   CMz_Mnt[iMarker] = 0.0;
      CoPx_Mnt[iMarker] = 0.0;         CoPy_Mnt[iMarker] = 0.0;  CoPz_Mnt[iMarker] = 0.0;
      CT_Mnt[iMarker] = 0.0;           CQ_Mnt[iMarker] = 0.0;    CMerit_Mnt[iMarker] = 0.0;
      CEff_Mnt[iMarker] = 0.0;
      
      for (iDim = 0; iDim < nDim; iDim++) ForceMomentum[iDim] = 0.0;
      MomentMomentum[0] = 0.0; MomentMomentum[1] = 0.0; MomentMomentum[2] = 0.0;
      MomentX_Force[0] = 0.0;  MomentX_Force[1] = 0.0;  MomentX_Force[2] = 0.0;
      MomentY_Force[0] = 0.0;  MomentY_Force[1] = 0.0;  MomentY_Force[2] = 0.0;
      MomentZ_Force[0] = 0.0;  MomentZ_Force[1] = 0.0;  MomentZ_Force[2] = 0.0;
      
      /*--- Loop over the vertices to compute the forces ---*/
      
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/
        
        if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {
          
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          Density   = node[iPoint]->GetDensity();
          
          /*--- Quadratic objective function for the near-field.
           This uses the infinity pressure regardless of Mach number. ---*/
          
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
          
          MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim]  = node[iPoint]->GetVelocity(iDim);
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
            MassFlow -= Normal[iDim]*Velocity[iDim]*Density;
          }
          
          /*--- Axisymmetric simulations ---*/
          
          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;
          
          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = MassFlow * Velocity[iDim] * factor * AxiFactor;
            ForceMomentum[iDim] += Force[iDim];
          }
          
          /*--- Moment with respect to the reference axis ---*/
          
          if (iDim == 3) {
            MomentMomentum[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLength;
            MomentX_Force[1]  += (-Force[1]*Coord[2]);
            MomentX_Force[2]  += (Force[2]*Coord[1]);
            
            MomentMomentum[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLength;
            MomentY_Force[2]  += (-Force[2]*Coord[0]);
            MomentY_Force[0]  += (Force[0]*Coord[2]);
          }
          MomentMomentum[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLength;
          MomentZ_Force[0]  += (-Force[0]*Coord[1]);
          MomentZ_Force[1]  += (Force[1]*Coord[0]);
          
        }
        
      }
      
      /*--- Project forces and store the non-dimensional coefficients ---*/
      
      if (Monitoring == YES) {
        
        if (nDim == 2) {
          CD_Mnt[iMarker]     =  ForceMomentum[0]*cos(Alpha) + ForceMomentum[1]*sin(Alpha);
          CL_Mnt[iMarker]     = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[1]*cos(Alpha);
          CEff_Mnt[iMarker]   = CL_Mnt[iMarker] / (CD_Mnt[iMarker]+EPS);
          CFx_Mnt[iMarker]    = ForceMomentum[0];
          CFy_Mnt[iMarker]    = ForceMomentum[1];
          CMz_Mnt[iMarker]    = MomentMomentum[2];
          CoPx_Mnt[iMarker]   = MomentZ_Force[1];
          CoPy_Mnt[iMarker]   = -MomentZ_Force[0];
          CT_Mnt[iMarker]     = -CFx_Mnt[iMarker];
          CQ_Mnt[iMarker]     = -CMz_Mnt[iMarker];
          CMerit_Mnt[iMarker] = CT_Mnt[iMarker] / (CQ_Mnt[iMarker] + EPS);
        }
        if (nDim == 3) {
          CD_Mnt[iMarker]         =  ForceMomentum[0]*cos(Alpha)*cos(Beta) + ForceMomentum[1]*sin(Beta) + ForceMomentum[2]*sin(Alpha)*cos(Beta);
          CL_Mnt[iMarker]         = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[2]*cos(Alpha);
          CSF_Mnt[iMarker]        = -ForceMomentum[0]*sin(Beta)*cos(Alpha) + ForceMomentum[1]*cos(Beta) - ForceMomentum[2]*sin(Beta)*sin(Alpha);
          CEff_Mnt[iMarker]       = CL_Mnt[iMarker] / (CD_Mnt[iMarker] + EPS);
          CFx_Mnt[iMarker]        = ForceMomentum[0];
          CFy_Mnt[iMarker]        = ForceMomentum[1];
          CFz_Mnt[iMarker]        = ForceMomentum[2];
          CMx_Mnt[iMarker]        = MomentMomentum[0];
          CMy_Mnt[iMarker]        = MomentMomentum[1];
          CMz_Mnt[iMarker]        = MomentMomentum[2];
          CoPx_Mnt[iMarker]       = -MomentY_Force[0];
          CoPz_Mnt[iMarker]       =  MomentY_Force[2];
          CT_Mnt[iMarker]         = -CFz_Mnt[iMarker];
          CQ_Mnt[iMarker]         = -CMz_Mnt[iMarker];
          CMerit_Mnt[iMarker]     = CT_Mnt[iMarker] / (CQ_Mnt[iMarker] + EPS);
        }
        
        AllBound_CD_Mnt           += CD_Mnt[iMarker];
        AllBound_CL_Mnt           += CL_Mnt[iMarker];
        AllBound_CSF_Mnt          += CSF_Mnt[iMarker];
        AllBound_CEff_Mnt          = AllBound_CL_Mnt / (AllBound_CD_Mnt + EPS);
        AllBound_CFx_Mnt          += CFx_Mnt[iMarker];
        AllBound_CFy_Mnt          += CFy_Mnt[iMarker];
        AllBound_CFz_Mnt          += CFz_Mnt[iMarker];
        AllBound_CMx_Mnt          += CMx_Mnt[iMarker];
        AllBound_CMy_Mnt          += CMy_Mnt[iMarker];
        AllBound_CMz_Mnt          += CMz_Mnt[iMarker];
        AllBound_CoPy_Mnt         += CoPy_Mnt[iMarker];
        AllBound_CoPz_Mnt         += CoPz_Mnt[iMarker];
        AllBound_CT_Mnt           += CT_Mnt[iMarker];
        AllBound_CQ_Mnt           += CQ_Mnt[iMarker];
        AllBound_CMerit_Mnt       += AllBound_CT_Mnt / (AllBound_CQ_Mnt + EPS);
        
        /*--- Compute the coefficients per surface ---*/
        
        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Mnt[iMarker_Monitoring]      += CL_Mnt[iMarker];
            Surface_CD_Mnt[iMarker_Monitoring]      += CD_Mnt[iMarker];
            Surface_CSF_Mnt[iMarker_Monitoring]     += CSF_Mnt[iMarker];
            Surface_CEff_Mnt[iMarker_Monitoring]     = CL_Mnt[iMarker] / (CD_Mnt[iMarker] + EPS);
            Surface_CFx_Mnt[iMarker_Monitoring]     += CFx_Mnt[iMarker];
            Surface_CFy_Mnt[iMarker_Monitoring]     += CFy_Mnt[iMarker];
            Surface_CFz_Mnt[iMarker_Monitoring]     += CFz_Mnt[iMarker];
            Surface_CMx_Mnt[iMarker_Monitoring]     += CMx_Mnt[iMarker];
            Surface_CMy_Mnt[iMarker_Monitoring]     += CMy_Mnt[iMarker];
            Surface_CMz_Mnt[iMarker_Monitoring]     += CMz_Mnt[iMarker];
          }
        }
        
      }
      
      
    }
  }
  
#ifdef HAVE_MPI
  
  /*--- Add AllBound information using all the nodes ---*/
  
  MyAllBound_CD_Mnt           = AllBound_CD_Mnt;           AllBound_CD_Mnt = 0.0;
  MyAllBound_CL_Mnt           = AllBound_CL_Mnt;           AllBound_CL_Mnt = 0.0;
  MyAllBound_CSF_Mnt          = AllBound_CSF_Mnt;          AllBound_CSF_Mnt = 0.0;
  MyAllBound_CFx_Mnt          = AllBound_CFx_Mnt;          AllBound_CFx_Mnt = 0.0;
  MyAllBound_CFy_Mnt          = AllBound_CFy_Mnt;          AllBound_CFy_Mnt = 0.0;
  MyAllBound_CFz_Mnt          = AllBound_CFz_Mnt;          AllBound_CFz_Mnt = 0.0;
  MyAllBound_CMx_Mnt          = AllBound_CMx_Mnt;          AllBound_CMx_Mnt = 0.0;
  MyAllBound_CMy_Mnt          = AllBound_CMy_Mnt;          AllBound_CMy_Mnt = 0.0;
  MyAllBound_CMz_Mnt          = AllBound_CMz_Mnt;          AllBound_CMz_Mnt = 0.0;
  MyAllBound_CoPx_Mnt         = AllBound_CoPx_Mnt;         AllBound_CoPx_Mnt = 0.0;
  MyAllBound_CoPy_Mnt         = AllBound_CoPy_Mnt;         AllBound_CoPy_Mnt = 0.0;
  MyAllBound_CoPz_Mnt         = AllBound_CoPz_Mnt;         AllBound_CoPz_Mnt = 0.0;
  MyAllBound_CT_Mnt           = AllBound_CT_Mnt;           AllBound_CT_Mnt = 0.0;
  MyAllBound_CQ_Mnt           = AllBound_CQ_Mnt;           AllBound_CQ_Mnt = 0.0;
  
  SU2_MPI::Allreduce(&MyAllBound_CD_Mnt, &AllBound_CD_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CL_Mnt, &AllBound_CL_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CSF_Mnt, &AllBound_CSF_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Mnt = AllBound_CL_Mnt / (AllBound_CD_Mnt + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CFx_Mnt, &AllBound_CFx_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFy_Mnt, &AllBound_CFy_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFz_Mnt, &AllBound_CFz_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMx_Mnt, &AllBound_CMx_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMy_Mnt, &AllBound_CMy_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMz_Mnt, &AllBound_CMz_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CoPx_Mnt, &AllBound_CoPx_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CoPy_Mnt, &AllBound_CoPy_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CoPz_Mnt, &AllBound_CoPz_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CT_Mnt, &AllBound_CT_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CQ_Mnt, &AllBound_CQ_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CMerit_Mnt = AllBound_CT_Mnt / (AllBound_CQ_Mnt + EPS);
  
  /*--- Add the forces on the surfaces using all the nodes ---*/
  
  MySurface_CL_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Mnt = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_CL_Mnt[iMarker_Monitoring]      = Surface_CL_Mnt[iMarker_Monitoring];
    MySurface_CD_Mnt[iMarker_Monitoring]      = Surface_CD_Mnt[iMarker_Monitoring];
    MySurface_CSF_Mnt[iMarker_Monitoring] = Surface_CSF_Mnt[iMarker_Monitoring];
    MySurface_CEff_Mnt[iMarker_Monitoring]       = Surface_CEff_Mnt[iMarker_Monitoring];
    MySurface_CFx_Mnt[iMarker_Monitoring]        = Surface_CFx_Mnt[iMarker_Monitoring];
    MySurface_CFy_Mnt[iMarker_Monitoring]        = Surface_CFy_Mnt[iMarker_Monitoring];
    MySurface_CFz_Mnt[iMarker_Monitoring]        = Surface_CFz_Mnt[iMarker_Monitoring];
    MySurface_CMx_Mnt[iMarker_Monitoring]        = Surface_CMx_Mnt[iMarker_Monitoring];
    MySurface_CMy_Mnt[iMarker_Monitoring]        = Surface_CMy_Mnt[iMarker_Monitoring];
    MySurface_CMz_Mnt[iMarker_Monitoring]        = Surface_CMz_Mnt[iMarker_Monitoring];
    
    Surface_CL_Mnt[iMarker_Monitoring]         = 0.0;
    Surface_CD_Mnt[iMarker_Monitoring]         = 0.0;
    Surface_CSF_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CEff_Mnt[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Mnt[iMarker_Monitoring]        = 0.0;
  }
  
  SU2_MPI::Allreduce(MySurface_CL_Mnt, Surface_CL_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CD_Mnt, Surface_CD_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CSF_Mnt, Surface_CSF_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
    Surface_CEff_Mnt[iMarker_Monitoring] = Surface_CL_Mnt[iMarker_Monitoring] / (Surface_CD_Mnt[iMarker_Monitoring] + EPS);
  SU2_MPI::Allreduce(MySurface_CFx_Mnt, Surface_CFx_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFy_Mnt, Surface_CFy_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFz_Mnt, Surface_CFz_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMx_Mnt, Surface_CMx_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMy_Mnt, Surface_CMy_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMz_Mnt, Surface_CMz_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  delete [] MySurface_CL_Mnt; delete [] MySurface_CD_Mnt; delete [] MySurface_CSF_Mnt;
  delete [] MySurface_CEff_Mnt;  delete [] MySurface_CFx_Mnt;   delete [] MySurface_CFy_Mnt;
  delete [] MySurface_CFz_Mnt;
  delete [] MySurface_CMx_Mnt;   delete [] MySurface_CMy_Mnt;  delete [] MySurface_CMz_Mnt;
  
#endif
  
  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/
  
  Total_CD            += AllBound_CD_Mnt;
  Total_CL            += AllBound_CL_Mnt;
  Total_CSF           += AllBound_CSF_Mnt;
  Total_CEff          = Total_CL / (Total_CD + EPS);
  Total_CFx           += AllBound_CFx_Mnt;
  Total_CFy           += AllBound_CFy_Mnt;
  Total_CFz           += AllBound_CFz_Mnt;
  Total_CMx           += AllBound_CMx_Mnt;
  Total_CMy           += AllBound_CMy_Mnt;
  Total_CMz           += AllBound_CMz_Mnt;
  Total_CoPx          += AllBound_CoPx_Mnt;
  Total_CoPy          += AllBound_CoPy_Mnt;
  Total_CoPz          += AllBound_CoPz_Mnt;
  Total_CT            += AllBound_CT_Mnt;
  Total_CQ            += AllBound_CQ_Mnt;
  Total_CMerit        = Total_CT / (Total_CQ + EPS);
  
  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]         += Surface_CL_Mnt[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]         += Surface_CD_Mnt[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring]        += Surface_CSF_Mnt[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       += Surface_CL_Mnt[iMarker_Monitoring] / (Surface_CD_Mnt[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        += Surface_CFx_Mnt[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        += Surface_CFy_Mnt[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        += Surface_CFz_Mnt[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        += Surface_CMx_Mnt[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        += Surface_CMy_Mnt[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        += Surface_CMz_Mnt[iMarker_Monitoring];
  }

}

void CEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {
  su2double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  bool adjoint = config->GetContinuous_Adjoint();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;
    
    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);
    
    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }
    
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
  
}

void CEulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {
  su2double *Residual, *Res_TruncError, Vol, Delta, Res, tmp_time, tmp_func;
  unsigned short iVar;
  unsigned long iPoint;

  /*--- Hard-coded classical RK4 coefficients. Will be added to config. ---*/
  su2double RK_FuncCoeff[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
  su2double RK_TimeCoeff[4] = {0.5, 0.5, 1.0, 1.0};

  bool adjoint = config->GetContinuous_Adjoint();

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;

    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);

    tmp_time = -1.0*RK_TimeCoeff[iRKStep]*Delta;
    tmp_func = -1.0*RK_FuncCoeff[iRKStep]*Delta;

    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        if (iRKStep < 3) {
          /* Base Solution Update */
          node[iPoint]->AddSolution(iVar, tmp_time*Res);

          /* New Solution Update */
          node[iPoint]->AddSolution_New(iVar, tmp_func*Res);
        } else {
          node[iPoint]->SetSolution(iVar, node[iPoint]->GetSolution_New(iVar) + tmp_func*Res);
        }

        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }

  }

  /*--- MPI solution ---*/

  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

void CEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  bool adjoint = config->GetContinuous_Adjoint();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    local_Residual = LinSysRes.GetBlock(iPoint);
    
    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = local_Residual[iVar] + local_Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }
    
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar, jVar;
  unsigned long iPoint, total_index, IterLinSol = 0;
  su2double Delta, *local_Res_TruncError, Vol;
  
  bool adjoint = config->GetContinuous_Adjoint();
  bool roe_turkel = config->GetKind_Upwind_Flow() == TURKEL;
  bool low_mach_prec = config->Low_Mach_Preconditioning();
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the residual ---*/
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    
    /*--- Read the volume ---*/
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    
    if (node[iPoint]->GetDelta_Time() != 0.0) {
      Delta = Vol / node[iPoint]->GetDelta_Time();
      if (roe_turkel || low_mach_prec) {
        SetPreconditioner(config, iPoint);
        for (iVar = 0; iVar < nVar; iVar ++ )
          for (jVar = 0; jVar < nVar; jVar ++ )
            LowMach_Precontioner[iVar][jVar] = Delta*LowMach_Precontioner[iVar][jVar];
        Jacobian.AddBlock(iPoint, iPoint, LowMach_Precontioner);
      }
      else {
        Jacobian.AddVal2Diag(iPoint, Delta);
      }
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }
  
  /*--- Solve or smooth the linear system ---*/
  
  CSysSolve system;
  IterLinSol = system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- The the number of iterations of the linear solver ---*/
  
  SetIterLinSolver(IterLinSol);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  if (!adjoint) {
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        node[iPoint]->AddSolution(iVar, config->GetRelaxation_Factor_Flow()*LinSysSol[iPoint*nVar+iVar]);
      }
    }
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CEulerSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iVar, iMarker;
  su2double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
  Partial_Gradient, Partial_Res, *Normal;
  
  /*--- Gradient primitive variables compressible (temp, vx, vy, vz, P, rho) ---*/

  PrimVar_Vertex = new su2double [nPrimVarGrad];
  PrimVar_i = new su2double [nPrimVarGrad];
  PrimVar_j = new su2double [nPrimVarGrad];
  
  /*--- Set Gradient_Primitive to zero ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);

  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
      PrimVar_j[iVar] = node[jPoint]->GetPrimitive(iVar);
    }
    
    Normal = geometry->edge[iEdge]->GetNormal();
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_Average =  0.5 * ( PrimVar_i[iVar] + PrimVar_j[iVar] );
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = PrimVar_Average*Normal[iDim];
        if (geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
        if (geometry->node[jPoint]->GetDomain())
          node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
      }
    }
  }

  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      if (geometry->node[iPoint]->GetDomain()) {
        
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          PrimVar_Vertex[iVar] = node[iPoint]->GetPrimitive(iVar);
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++) {
            Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
            node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
          }
      }
    }
  }

  /*--- Update gradient value ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Gradient = node[iPoint]->GetGradient_Primitive(iVar, iDim) / (geometry->node[iPoint]->GetVolume());
        node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);
      }
    }
  }

  delete [] PrimVar_Vertex;
  delete [] PrimVar_i;
  delete [] PrimVar_j;

  Set_MPI_Primitive_Gradient(geometry, config);

}

void CEulerSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
  bool singular;
  
  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Set the value of the singular ---*/
    singular = false;
    
    /*--- Get coordinates ---*/
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    
    /*--- Get primitives from CVariable ---*/
    
    PrimVar_i = node[iPoint]->GetPrimitive();
    
    /*--- Inizialization of variables ---*/
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Cvector[iVar][iDim] = 0.0;
    
    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
    
    AD::StartPreacc();
    AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    AD::SetPreaccIn(Coord_i, nDim);
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      
      PrimVar_j = node[jPoint]->GetPrimitive();
      
      AD::SetPreaccIn(Coord_j, nDim);
      AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);

      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      
      if (weight != 0.0) {
        
        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
        
        if (nDim == 3) {
          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
        }
        
        /*--- Entries of c:= transpose(A)*b ---*/
        
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
        
      }
      
    }
    
    /*--- Entries of upper triangular matrix R ---*/
    
    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
    
    if (nDim == 3) {
      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
    }
    
    /*--- Compute determinant ---*/
    
    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);
    
    /*--- Detect singular matrices ---*/
    
    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }
    
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    
    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }
    
    /*--- Computation of the gradient: S*c ---*/
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        product = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
        }
        
        node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
      }
    }
    
    AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
    AD::EndPreacc();
  }
  
  Set_MPI_Primitive_Gradient(geometry, config);
  
}

void CEulerSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  su2double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j,
  *Primitive, *Primitive_i, *Primitive_j, *LocalMinPrimitive, *LocalMaxPrimitive,
  *GlobalMinPrimitive, *GlobalMaxPrimitive,
  dave, LimK, eps2, eps1, dm, dp, du, y, limiter;
  
  dave = config->GetRefElemLength();
  LimK = config->GetVenkat_LimiterCoeff();

  if (config->GetKind_SlopeLimit_Flow() == NO_LIMITER) {
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        node[iPoint]->SetLimiter_Primitive(iVar, 1.0);
      }
    }
    
  }
  
  else {
    
    /*--- Initialize solution max and solution min and the limiter in the entire domain --*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        node[iPoint]->SetSolution_Max(iVar, -EPS);
        node[iPoint]->SetSolution_Min(iVar, EPS);
        node[iPoint]->SetLimiter_Primitive(iVar, 2.0);
      }
    }
    
    /*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      
      /*--- Get the primitive variables ---*/
      
      Primitive_i = node[iPoint]->GetPrimitive();
      Primitive_j = node[jPoint]->GetPrimitive();
      
      /*--- Compute the maximum, and minimum values for nodes i & j ---*/
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        du = (Primitive_j[iVar] - Primitive_i[iVar]);
        node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
        node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
        node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
        node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
      }
      
    }
    
  }
  
  /*--- Barth-Jespersen limiter with Venkatakrishnan modification ---*/
  
  if (config->GetKind_SlopeLimit_Flow() == BARTH_JESPERSEN) {
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        
        AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
          else dp = node[iPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }
        
        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
        }
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
          else dp = node[jPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }
        
        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }
        
      }
      
      AD::EndPreacc();

    }
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        y =  node[iPoint]->GetLimiter_Primitive(iVar);
        limiter = (y*y + 2.0*y) / (y*y + y + 2.0);
        node[iPoint]->SetLimiter_Primitive(iVar, limiter);
      }
    }
    
  }
  
  /*--- Venkatakrishnan limiter ---*/
  
  if ((config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN) ||
      (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG)) {
    
    /*--- Allocate memory for the max and min primitive value --*/
    
    LocalMinPrimitive = new su2double [nPrimVarGrad]; GlobalMinPrimitive = new su2double [nPrimVarGrad];
    LocalMaxPrimitive = new su2double [nPrimVarGrad]; GlobalMaxPrimitive = new su2double [nPrimVarGrad];
    
    /*--- Compute the max value and min value of the solution ---*/
    
    Primitive = node[0]->GetPrimitive();
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      LocalMinPrimitive[iVar] = Primitive[iVar];
      LocalMaxPrimitive[iVar] = Primitive[iVar];
    }
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Get the primitive variables ---*/
      
      Primitive = node[iPoint]->GetPrimitive();
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        LocalMinPrimitive[iVar] = min (LocalMinPrimitive[iVar], Primitive[iVar]);
        LocalMaxPrimitive[iVar] = max (LocalMaxPrimitive[iVar], Primitive[iVar]);
      }
      
    }
    
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(LocalMinPrimitive, GlobalMinPrimitive, nPrimVarGrad, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(LocalMaxPrimitive, GlobalMaxPrimitive, nPrimVarGrad, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      GlobalMinPrimitive[iVar] = LocalMinPrimitive[iVar];
      GlobalMaxPrimitive[iVar] = LocalMaxPrimitive[iVar];
    }
#endif
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
      AD::SetPreaccIn(eps2);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        
        if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG) {
          eps1 = LimK * (GlobalMaxPrimitive[iVar] - GlobalMinPrimitive[iVar]);
          eps2 = eps1*eps1;
        }
        else {
          eps1 = LimK*dave;
          eps2 = eps1*eps1*eps1;
        }

        AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
        else dp = node[iPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
        }
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }
        
      }

      AD::EndPreacc();
      
    }
    
    delete [] LocalMinPrimitive; delete [] GlobalMinPrimitive;
    delete [] LocalMaxPrimitive; delete [] GlobalMaxPrimitive;
    
  }

  /*--- Limiter MPI ---*/
  
  Set_MPI_Primitive_Limiter(geometry, config);

}

void CEulerSolver::SetPreconditioner(CConfig *config, unsigned long iPoint) {
  unsigned short iDim, jDim, iVar, jVar;
  su2double local_Mach, rho, enthalpy, soundspeed, sq_vel;
  su2double *U_i = NULL;
  su2double Beta_max = config->GetmaxTurkelBeta();
  su2double Mach_infty2, Mach_lim2, aux, parameter;
  
  /*--- Variables to calculate the preconditioner parameter Beta ---*/
  local_Mach = sqrt(node[iPoint]->GetVelocity2())/node[iPoint]->GetSoundSpeed();

  /*--- Weiss and Smith Preconditioning---*/
  Mach_infty2 = pow(config->GetMach(),2.0);
  Mach_lim2 = pow(0.00001,2.0);
  aux = max(pow(local_Mach,2.0),Mach_lim2);
  parameter = min(1.0, max(aux,Beta_max*Mach_infty2));
  
  U_i = node[iPoint]->GetSolution();
  
  rho = U_i[0];
  enthalpy = node[iPoint]->GetEnthalpy();
  soundspeed = node[iPoint]->GetSoundSpeed();
  sq_vel = node[iPoint]->GetVelocity2();
  
  /*---Calculating the inverse of the preconditioning matrix that multiplies the time derivative  */
  LowMach_Precontioner[0][0] = 0.5*sq_vel;
  LowMach_Precontioner[0][nVar-1] = 1.0;
  for (iDim = 0; iDim < nDim; iDim ++)
    LowMach_Precontioner[0][1+iDim] = -1.0*U_i[iDim+1]/rho;
  
  for (iDim = 0; iDim < nDim; iDim ++) {
    LowMach_Precontioner[iDim+1][0] = 0.5*sq_vel*U_i[iDim+1]/rho;
    LowMach_Precontioner[iDim+1][nVar-1] = U_i[iDim+1]/rho;
    for (jDim = 0; jDim < nDim; jDim ++) {
      LowMach_Precontioner[iDim+1][1+jDim] = -1.0*U_i[jDim+1]/rho*U_i[iDim+1]/rho;
    }
  }
  
  LowMach_Precontioner[nVar-1][0] = 0.5*sq_vel*enthalpy;
  LowMach_Precontioner[nVar-1][nVar-1] = enthalpy;
  for (iDim = 0; iDim < nDim; iDim ++)
    LowMach_Precontioner[nVar-1][1+iDim] = -1.0*U_i[iDim+1]/rho*enthalpy;
  
  
  for (iVar = 0; iVar < nVar; iVar ++ ) {
    for (jVar = 0; jVar < nVar; jVar ++ ) {
      LowMach_Precontioner[iVar][jVar] = (parameter - 1.0) * ((Gamma-1.0)/(soundspeed*soundspeed))*LowMach_Precontioner[iVar][jVar];
      if (iVar == jVar)
        LowMach_Precontioner[iVar][iVar] += 1.0;
    }
  }
  
}

void CEulerSolver::GetPower_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {
  
  unsigned short iDim, iMarker, jMarker;
  unsigned long iVertex, iPoint;
  su2double  *V_inlet = NULL, *V_outlet = NULL, Pressure, Temperature, Velocity[3], Vn,
  Velocity2, Density, Area, SoundSpeed, TotalPressure, Vel_Infty2, RamDrag,
  TotalTemperature, VelocityJet,
  Vel_Infty, MaxPressure, MinPressure, MFR, InfVel2;
  unsigned short iMarker_Inlet, iMarker_Outlet, nMarker_Inlet, nMarker_Outlet;
  string Inlet_TagBound, Outlet_TagBound;
  su2double DeltaPress = 0.0, DeltaTemp = 0.0, TotalPressRatio = 0.0, TotalTempRatio = 0.0, StaticPressRatio = 0.0, StaticTempRatio = 0.0,
  NetThrust = 0.0, GrossThrust = 0.0, Power = 0.0, MassFlow = 0.0, Mach = 0.0, Force = 0.0;
  bool ReverseFlow, Engine = false, Pair = true;
  
  su2double Gas_Constant                         = config->GetGas_ConstantND();
  su2double Cp                                               = Gas_Constant*Gamma / (Gamma-1.0);
  su2double Alpha                                         = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta                                           = config->GetAoS()*PI_NUMBER/180.0;
  bool write_heads = ((((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0) && (config->GetExtIter()!= 0))
                      || (config->GetExtIter() == 1));
  bool Evaluate_BC = ((((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0))
                      || (config->GetExtIter() == 1) || (config->GetDiscrete_Adjoint()));
  
  if ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0)) Engine = true;
  if ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0)) Engine = false;
  if ((config->GetnMarker_EngineInflow()) != (config->GetnMarker_EngineExhaust())) Pair = false;
  
  
  if (Engine) { nMarker_Inlet  = config->GetnMarker_EngineInflow(); nMarker_Outlet = config->GetnMarker_EngineExhaust(); }
  else  { nMarker_Inlet   = config->GetnMarker_ActDiskInlet(); nMarker_Outlet  = config->GetnMarker_ActDiskOutlet(); }
  
  /*--- Evaluate the MPI for the actuator disk IO ---*/
  
  if (Evaluate_BC) {
    
    /*--- Allocate memory ---*/
    
    su2double *Inlet_MassFlow         = new su2double [config->GetnMarker_All()];
    su2double *Inlet_ReverseMassFlow  = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Pressure         = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Mach             = new su2double [config->GetnMarker_All()];
    su2double *Inlet_MaxPressure      = new su2double [config->GetnMarker_All()];
    su2double *Inlet_MinPressure      = new su2double [config->GetnMarker_All()];
    su2double *Inlet_TotalPressure    = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Temperature      = new su2double [config->GetnMarker_All()];
    su2double *Inlet_TotalTemperature   = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Area             = new su2double [config->GetnMarker_All()];
    su2double *Inlet_RamDrag          = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Force            = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Power            = new su2double [config->GetnMarker_All()];
    su2double *Inlet_XCG                        = new su2double [config->GetnMarker_All()];
    su2double *Inlet_YCG                        = new su2double [config->GetnMarker_All()];
    su2double *Inlet_ZCG                        = new su2double [config->GetnMarker_All()];
    
    su2double *Outlet_MassFlow           = new su2double [config->GetnMarker_All()];
    su2double *Outlet_Pressure           = new su2double [config->GetnMarker_All()];
    su2double *Outlet_TotalPressure    = new su2double [config->GetnMarker_All()];
    su2double *Outlet_Temperature        = new su2double [config->GetnMarker_All()];
    su2double *Outlet_TotalTemperature = new su2double [config->GetnMarker_All()];
    su2double *Outlet_Area                   = new su2double [config->GetnMarker_All()];
    su2double *Outlet_GrossThrust        = new su2double [config->GetnMarker_All()];
    su2double *Outlet_Force                    = new su2double [config->GetnMarker_All()];
    su2double *Outlet_Power                    = new su2double [config->GetnMarker_All()];
    
    /*--- Comute MassFlow, average temp, press, etc. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      Inlet_MassFlow[iMarker] = 0.0;      Inlet_ReverseMassFlow[iMarker] = 0.0;  Inlet_MinPressure[iMarker] = 0.0;
      Inlet_Pressure[iMarker] = 0.0;      Inlet_Mach[iMarker] = 0.0;             Inlet_Temperature[iMarker] = 0.0;
      Inlet_MinPressure[iMarker] = +1E10; Inlet_MaxPressure[iMarker] = -1E10;    Inlet_Power[iMarker] = 0.0;
      Inlet_TotalPressure[iMarker] = 0.0; Inlet_TotalTemperature[iMarker] = 0.0;
      Inlet_Area[iMarker] = 0.0;
      Inlet_RamDrag[iMarker] = 0.0; Inlet_Force[iMarker]    = 0.0;
      Inlet_XCG[iMarker] = 0.0;  Inlet_YCG[iMarker]    = 0.0; Inlet_ZCG[iMarker]    = 0.0;
      
      Outlet_MassFlow[iMarker] = 0.0;
      Outlet_Pressure[iMarker] = 0.0; Outlet_Temperature[iMarker] = 0.0;
      Outlet_TotalPressure[iMarker] = 0.0; Outlet_TotalTemperature[iMarker] = 0.0;
      Outlet_Area[iMarker] = 0.0;
      Outlet_GrossThrust[iMarker]    = 0.0; Outlet_Force[iMarker]    = 0.0; Outlet_Power[iMarker]    = 0.0;
      
      MinPressure = +1E10; MaxPressure = -1E10;
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (geometry->node[iPoint]->GetDomain()) {
            
            V_inlet = node[iPoint]->GetPrimitive();
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
            
            Temperature     = V_inlet[0];
            Pressure                = V_inlet[nDim+1];
            
            Density                         = V_inlet[nDim+2];
            SoundSpeed  = sqrt(Gamma*Pressure/Density);
            
            Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vel_Infty2 =0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Area += Vector[iDim]*Vector[iDim];
              Velocity[iDim] = V_inlet[iDim+1];
              Velocity2 += Velocity[iDim]*Velocity[iDim];
              Vel_Infty2 += GetVelocity_Inf(iDim)*GetVelocity_Inf(iDim);
              MassFlow -= Vector[iDim]*Velocity[iDim]*Density;
            }
            
            Vn = 0.0; ReverseFlow = false;
            for (iDim = 0; iDim < nDim; iDim++) {  Vn -= Velocity[iDim]*Vector[iDim]/Area; }
            if (Vn < 0.0) { ReverseFlow = true; }
            
            Vel_Infty                   = sqrt (Vel_Infty2);
            Area                                = sqrt (Area);
            Mach                                = sqrt(Velocity2)/SoundSpeed;
            TotalPressure           = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma    / (Gamma - 1.0));
            TotalTemperature  = Temperature * (1.0 + Mach * Mach * 0.5 * (Gamma - 1.0));
            MinPressure               = min(MinPressure, TotalPressure);
            MaxPressure               = max(MaxPressure, TotalPressure);
            
            RamDrag     = MassFlow * Vel_Infty;
            
            Inlet_MassFlow[iMarker]         += MassFlow;
            Inlet_Pressure[iMarker]         += Pressure*MassFlow;
            Inlet_Mach[iMarker]             += Mach*MassFlow;
            Inlet_MinPressure[iMarker]       = min (MinPressure, Inlet_MinPressure[iMarker]);
            Inlet_MaxPressure[iMarker]       = max(MaxPressure, Inlet_MaxPressure[iMarker]);
            Inlet_TotalPressure[iMarker]    += TotalPressure*MassFlow;
            Inlet_Temperature[iMarker]      += Temperature*MassFlow;
            Inlet_TotalTemperature[iMarker] += TotalTemperature*MassFlow;
            Inlet_Area[iMarker]             += Area;
            Inlet_RamDrag[iMarker]          += RamDrag;
            Inlet_Power[iMarker] += MassFlow*Cp*TotalTemperature;
            if (ReverseFlow) Inlet_ReverseMassFlow[iMarker]  += MassFlow;
            
            su2double Inlet_ForceX =  -(Pressure - Pressure_Inf)*Vector[0] + MassFlow*Velocity[0];
            su2double Inlet_ForceY = -(Pressure - Pressure_Inf)*Vector[1] + MassFlow*Velocity[1];
            su2double Inlet_ForceZ = 0.0;
            if (nDim == 3) Inlet_ForceZ = -(Pressure - Pressure_Inf)*Vector[2] + MassFlow*Velocity[2];
            Inlet_Force[iMarker] +=  Inlet_ForceX*cos(Alpha)*cos(Beta) + Inlet_ForceY*sin(Beta) +Inlet_ForceZ*sin(Alpha)*cos(Beta);
            
            Inlet_XCG[iMarker]                        += geometry->node[iPoint]->GetCoord(0)*Area;
            Inlet_YCG[iMarker]                        += geometry->node[iPoint]->GetCoord(1)*Area;
            if (nDim == 3) Inlet_ZCG[iMarker] += geometry->node[iPoint]->GetCoord(2)*Area;
            
          }
        }
        
      }
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (geometry->node[iPoint]->GetDomain()) {
            
            V_outlet = node[iPoint]->GetPrimitive();
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
            
            Temperature         = V_outlet[0];
            Pressure                    = V_outlet[nDim+1];
            
            Density                         = V_outlet[nDim+2];
            SoundSpeed  = sqrt(Gamma*Pressure/Density);
            
            Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vel_Infty2 = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Area += Vector[iDim]*Vector[iDim];
              Velocity[iDim] = V_outlet[iDim+1];
              Velocity2 += Velocity[iDim]*Velocity[iDim];
              Vel_Infty2 += GetVelocity_Inf(iDim)*GetVelocity_Inf(iDim);
              MassFlow += Vector[iDim]*Velocity[iDim]*Density;
            }
            
            Vel_Infty                   = sqrt (Vel_Infty2);
            Area                                            = sqrt (Area);
            Mach                                = sqrt(Velocity2)/SoundSpeed;
            TotalPressure           = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma    / (Gamma - 1.0));
            TotalTemperature    = Temperature * (1.0 + Mach * Mach * 0.5 * (Gamma - 1.0));
            VelocityJet         = sqrt(Velocity2) ;

            GrossThrust     = MassFlow * VelocityJet;

            Outlet_MassFlow[iMarker]                            += MassFlow;
            Outlet_Pressure[iMarker]                                += Pressure*MassFlow;
            Outlet_TotalPressure[iMarker]               += TotalPressure*MassFlow;
            Outlet_Temperature[iMarker]                     += Temperature*MassFlow;
            Outlet_TotalTemperature[iMarker]    += TotalTemperature*MassFlow;
            Outlet_Area[iMarker]                                            += Area;
            Outlet_GrossThrust[iMarker]             += GrossThrust;
            Outlet_Power[iMarker] += MassFlow*Cp*TotalTemperature;
            
            su2double Outlet_ForceX = -(Pressure - Pressure_Inf)*Vector[0] -MassFlow*Velocity[0];
            su2double Outlet_ForceY =  -(Pressure - Pressure_Inf)*Vector[1] -MassFlow*Velocity[1];
            su2double Outlet_ForceZ = 0.0;
            if (nDim == 3) Outlet_ForceZ = -(Pressure - Pressure_Inf)*Vector[2] -MassFlow*Velocity[2];
            
            if (nDim == 2) Outlet_Force[iMarker] +=  Outlet_ForceX*cos(Alpha) + Outlet_ForceY*sin(Alpha);
            if (nDim == 3) Outlet_Force[iMarker] +=  Outlet_ForceX*cos(Alpha)*cos(Beta) + Outlet_ForceY*sin(Beta) + Outlet_ForceZ*sin(Alpha)*cos(Beta);
            
          }
        }
        
      }
      
    }
    
    /*--- Copy to the appropriate structure ---*/
    
    su2double *Inlet_MassFlow_Local             = new su2double [nMarker_Inlet];
    su2double *Inlet_ReverseMassFlow_Local              = new su2double [nMarker_Inlet];
    su2double *Inlet_Temperature_Local          = new su2double [nMarker_Inlet];
    su2double *Inlet_TotalTemperature_Local = new su2double [nMarker_Inlet];
    su2double *Inlet_Pressure_Local             = new su2double [nMarker_Inlet];
    su2double *Inlet_Mach_Local             = new su2double [nMarker_Inlet];
    su2double *Inlet_MinPressure_Local      = new su2double [nMarker_Inlet];
    su2double *Inlet_MaxPressure_Local      = new su2double [nMarker_Inlet];
    su2double *Inlet_Power_Local        = new su2double [nMarker_Inlet];
    su2double *Inlet_TotalPressure_Local        = new su2double [nMarker_Inlet];
    su2double *Inlet_RamDrag_Local              = new su2double [nMarker_Inlet];
    su2double *Inlet_Force_Local                            = new su2double [nMarker_Inlet];
    su2double *Inlet_Area_Local              = new su2double [nMarker_Inlet];
    su2double *Inlet_XCG_Local                      = new su2double [nMarker_Inlet];
    su2double *Inlet_YCG_Local                      = new su2double [nMarker_Inlet];
    su2double *Inlet_ZCG_Local                      = new su2double [nMarker_Inlet];
    
    su2double *Inlet_MassFlow_Total                         = new su2double [nMarker_Inlet];
    su2double *Inlet_ReverseMassFlow_Total              = new su2double [nMarker_Inlet];
    su2double *Inlet_Pressure_Total                             = new su2double [nMarker_Inlet];
    su2double *Inlet_Mach_Total                             = new su2double [nMarker_Inlet];
    su2double *Inlet_MinPressure_Total            = new su2double [nMarker_Inlet];
    su2double *Inlet_MaxPressure_Total         = new su2double [nMarker_Inlet];
    su2double *Inlet_Power_Total                            = new su2double [nMarker_Inlet];
    su2double *Inlet_TotalPressure_Total            = new su2double [nMarker_Inlet];
    su2double *Inlet_Temperature_Total                  = new su2double [nMarker_Inlet];
    su2double *Inlet_TotalTemperature_Total     = new su2double [nMarker_Inlet];
    su2double *Inlet_RamDrag_Total                          = new su2double [nMarker_Inlet];
    su2double *Inlet_Force_Total                                = new su2double [nMarker_Inlet];
    su2double *Inlet_Area_Total                                         = new su2double [nMarker_Inlet];
    su2double *Inlet_XCG_Total                                          = new su2double [nMarker_Inlet];
    su2double *Inlet_YCG_Total                                          = new su2double [nMarker_Inlet];
    su2double *Inlet_ZCG_Total                                          = new su2double [nMarker_Inlet];
    
    for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
      Inlet_MassFlow_Local[iMarker_Inlet]                       = 0.0;
      Inlet_ReverseMassFlow_Local[iMarker_Inlet]            = 0.0;
      Inlet_Pressure_Local[iMarker_Inlet]                           = 0.0;
      Inlet_Mach_Local[iMarker_Inlet]                           = 0.0;
      Inlet_MinPressure_Local[iMarker_Inlet]            = 0.0;
      Inlet_MaxPressure_Local[iMarker_Inlet]            = 0.0;
      Inlet_TotalPressure_Local[iMarker_Inlet]          = 0.0;
      Inlet_Temperature_Local[iMarker_Inlet]                    = 0.0;
      Inlet_TotalTemperature_Local[iMarker_Inlet] = 0.0;
      Inlet_RamDrag_Local[iMarker_Inlet]                    = 0.0;
      Inlet_Force_Local[iMarker_Inlet]                              = 0.0;
      Inlet_Power_Local[iMarker_Inlet]                          = 0.0;
      Inlet_Area_Local[iMarker_Inlet]                                           = 0.0;
      Inlet_XCG_Local[iMarker_Inlet]                                            = 0.0;
      Inlet_YCG_Local[iMarker_Inlet]                                            = 0.0;
      Inlet_ZCG_Local[iMarker_Inlet]                                            = 0.0;
      
      Inlet_MassFlow_Total[iMarker_Inlet]                           = 0.0;
      Inlet_ReverseMassFlow_Total[iMarker_Inlet]                = 0.0;
      Inlet_Pressure_Total[iMarker_Inlet]                               = 0.0;
      Inlet_Mach_Total[iMarker_Inlet]                               = 0.0;
      Inlet_MinPressure_Total[iMarker_Inlet]                = 0.0;
      Inlet_MaxPressure_Total[iMarker_Inlet]                = 0.0;
      Inlet_TotalPressure_Total[iMarker_Inlet]              = 0.0;
      Inlet_Temperature_Total[iMarker_Inlet]                    = 0.0;
      Inlet_TotalTemperature_Total[iMarker_Inlet]   = 0.0;
      Inlet_RamDrag_Total[iMarker_Inlet]                        = 0.0;
      Inlet_Force_Total[iMarker_Inlet]                                  = 0.0;
      Inlet_Power_Total[iMarker_Inlet]              = 0.0;
      Inlet_Area_Total[iMarker_Inlet]                                           = 0.0;
      Inlet_XCG_Total[iMarker_Inlet]                                            = 0.0;
      Inlet_YCG_Total[iMarker_Inlet]                                            = 0.0;
      Inlet_ZCG_Total[iMarker_Inlet]                                            = 0.0;
    }
    
    su2double *Outlet_MassFlow_Local                            = new su2double [nMarker_Outlet];
    su2double *Outlet_Pressure_Local                                = new su2double [nMarker_Outlet];
    su2double *Outlet_TotalPressure_Local           = new su2double [nMarker_Outlet];
    su2double *Outlet_Temperature_Local                     = new su2double [nMarker_Outlet];
    su2double *Outlet_TotalTemperature_Local    = new su2double [nMarker_Outlet];
    su2double *Outlet_GrossThrust_Local             = new su2double [nMarker_Outlet];
    su2double *Outlet_Force_Local               = new su2double [nMarker_Outlet];
    su2double *Outlet_Power_Local               = new su2double [nMarker_Outlet];
    su2double *Outlet_Area_Local                                            = new su2double [nMarker_Outlet];
    
    su2double *Outlet_MassFlow_Total                        = new su2double [nMarker_Outlet];
    su2double *Outlet_Pressure_Total                            = new su2double [nMarker_Outlet];
    su2double *Outlet_TotalPressure_Total           = new su2double [nMarker_Outlet];
    su2double *Outlet_Temperature_Total                     = new su2double [nMarker_Outlet];
    su2double *Outlet_TotalTemperature_Total    = new su2double [nMarker_Outlet];
    su2double *Outlet_GrossThrust_Total             = new su2double [nMarker_Outlet];
    su2double *Outlet_Force_Total                             = new su2double [nMarker_Outlet];
    su2double *Outlet_Power_Total                             = new su2double [nMarker_Outlet];
    su2double *Outlet_Area_Total                                            = new su2double [nMarker_Outlet];
    
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      Outlet_MassFlow_Local[iMarker_Outlet]                     = 0.0;
      Outlet_Pressure_Local[iMarker_Outlet]                         = 0.0;
      Outlet_TotalPressure_Local[iMarker_Outlet]            = 0.0;
      Outlet_Temperature_Local[iMarker_Outlet]                  = 0.0;
      Outlet_TotalTemperature_Local[iMarker_Outlet] = 0.0;
      Outlet_GrossThrust_Local[iMarker_Outlet]              = 0.0;
      Outlet_Force_Local[iMarker_Outlet]                        = 0.0;
      Outlet_Power_Local[iMarker_Outlet]                        = 0.0;
      Outlet_Area_Local[iMarker_Outlet]                                         = 0.0;
      
      Outlet_MassFlow_Total[iMarker_Outlet]                             = 0.0;
      Outlet_Pressure_Total[iMarker_Outlet]                             = 0.0;
      Outlet_TotalPressure_Total[iMarker_Outlet]                = 0.0;
      Outlet_Temperature_Total[iMarker_Outlet]                  = 0.0;
      Outlet_TotalTemperature_Total[iMarker_Outlet] = 0.0;
      Outlet_GrossThrust_Total[iMarker_Outlet]              = 0.0;
      Outlet_Force_Total[iMarker_Outlet]                        = 0.0;
      Outlet_Power_Total[iMarker_Outlet]                        = 0.0;
      Outlet_Area_Total[iMarker_Outlet]                                         = 0.0;
    }
    
    /*--- Copy the values to the local array for MPI ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||  (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW)) {
        for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
          
          if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) Inlet_TagBound = config->GetMarker_ActDiskInlet_TagBound(iMarker_Inlet);
          if (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW) Inlet_TagBound = config->GetMarker_EngineInflow_TagBound(iMarker_Inlet);
          
          if (config->GetMarker_All_TagBound(iMarker) == Inlet_TagBound) {
            Inlet_MassFlow_Local[iMarker_Inlet]                             += Inlet_MassFlow[iMarker];
            Inlet_ReverseMassFlow_Local[iMarker_Inlet]              += Inlet_ReverseMassFlow[iMarker];
            Inlet_Pressure_Local[iMarker_Inlet]                                 += Inlet_Pressure[iMarker];
            Inlet_Mach_Local[iMarker_Inlet]                                     += Inlet_Mach[iMarker];
            Inlet_MinPressure_Local[iMarker_Inlet]              += Inlet_MinPressure[iMarker];
            Inlet_MaxPressure_Local[iMarker_Inlet]          += Inlet_MaxPressure[iMarker];
            Inlet_TotalPressure_Local[iMarker_Inlet]            += Inlet_TotalPressure[iMarker];
            Inlet_Temperature_Local[iMarker_Inlet]                  += Inlet_Temperature[iMarker];
            Inlet_TotalTemperature_Local[iMarker_Inlet] += Inlet_TotalTemperature[iMarker];
            Inlet_RamDrag_Local[iMarker_Inlet]                          += Inlet_RamDrag[iMarker];
            Inlet_Force_Local[iMarker_Inlet]                                    += Inlet_Force[iMarker];
            Inlet_Power_Local[iMarker_Inlet]                          += Inlet_Power[iMarker];
            Inlet_Area_Local[iMarker_Inlet]                                                 += Inlet_Area[iMarker];
            Inlet_XCG_Local[iMarker_Inlet]                                              += Inlet_XCG[iMarker];
            Inlet_YCG_Local[iMarker_Inlet]                                              += Inlet_YCG[iMarker];
            if (nDim == 3) Inlet_ZCG_Local[iMarker_Inlet]                                               += Inlet_ZCG[iMarker];
          }
          
        }
      }
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) || (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST)) {
        for (iMarker_Outlet= 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
          
          if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) Outlet_TagBound = config->GetMarker_ActDiskOutlet_TagBound(iMarker_Outlet);
          if (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST) Outlet_TagBound = config->GetMarker_EngineExhaust_TagBound(iMarker_Outlet);
          
          if (config->GetMarker_All_TagBound(iMarker) == Outlet_TagBound) {
            Outlet_MassFlow_Local[iMarker_Outlet]                           += Outlet_MassFlow[iMarker];
            Outlet_Pressure_Local[iMarker_Outlet]                               += Outlet_Pressure[iMarker];
            Outlet_TotalPressure_Local[iMarker_Outlet]          += Outlet_TotalPressure[iMarker];
            Outlet_Temperature_Local[iMarker_Outlet]                    += Outlet_Temperature[iMarker];
            Outlet_TotalTemperature_Local[iMarker_Outlet]   += Outlet_TotalTemperature[iMarker];
            Outlet_GrossThrust_Local[iMarker_Outlet]                    += Outlet_GrossThrust[iMarker];
            Outlet_Force_Local[iMarker_Outlet]                              += Outlet_Force[iMarker];
            Outlet_Power_Local[iMarker_Outlet]                              += Outlet_Power[iMarker];
            Outlet_Area_Local[iMarker_Outlet]                                               += Outlet_Area[iMarker];
          }
          
        }
      }
      
    }
    
    /*--- Correct the min max values for the MPI ---*/
    
    bool ActDisk  = false;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW)) { ActDisk  = true; break; }
    }
    
    if (!ActDisk) {
      for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
        Inlet_MinPressure_Local[iMarker_Inlet]          = 1E10;
        Inlet_MaxPressure_Local[iMarker_Inlet]      = -1E10;
      }
    }
    
    /*--- All the ranks to compute the total value ---*/
    
#ifdef HAVE_MPI
    
    SU2_MPI::Allreduce(Inlet_MassFlow_Local, Inlet_MassFlow_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_ReverseMassFlow_Local, Inlet_ReverseMassFlow_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Pressure_Local, Inlet_Pressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Mach_Local, Inlet_Mach_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_MinPressure_Local, Inlet_MinPressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_MaxPressure_Local, Inlet_MaxPressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_TotalPressure_Local, Inlet_TotalPressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Temperature_Local, Inlet_Temperature_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_TotalTemperature_Local, Inlet_TotalTemperature_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_RamDrag_Local, Inlet_RamDrag_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Force_Local, Inlet_Force_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Power_Local, Inlet_Power_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Area_Local, Inlet_Area_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_XCG_Local, Inlet_XCG_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_YCG_Local, Inlet_YCG_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (nDim == 3) SU2_MPI::Allreduce(Inlet_ZCG_Local, Inlet_ZCG_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    SU2_MPI::Allreduce(Outlet_MassFlow_Local, Outlet_MassFlow_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Pressure_Local, Outlet_Pressure_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_TotalPressure_Local, Outlet_TotalPressure_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Temperature_Local, Outlet_Temperature_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_TotalTemperature_Local, Outlet_TotalTemperature_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_GrossThrust_Local, Outlet_GrossThrust_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Force_Local, Outlet_Force_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Power_Local, Outlet_Power_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Area_Local, Outlet_Area_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
#else
    
    for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
      Inlet_MassFlow_Total[iMarker_Inlet]                            = Inlet_MassFlow_Local[iMarker_Inlet];
      Inlet_ReverseMassFlow_Total[iMarker_Inlet] = Inlet_ReverseMassFlow_Local[iMarker_Inlet];
      Inlet_Pressure_Total[iMarker_Inlet]                              = Inlet_Pressure_Local[iMarker_Inlet];
      Inlet_Mach_Total[iMarker_Inlet]                                        = Inlet_Mach_Local[iMarker_Inlet];
      Inlet_MinPressure_Total[iMarker_Inlet]                 = Inlet_MinPressure_Local[iMarker_Inlet];
      Inlet_MaxPressure_Total[iMarker_Inlet]             = Inlet_MaxPressure_Local[iMarker_Inlet];
      Inlet_TotalPressure_Total[iMarker_Inlet]               = Inlet_TotalPressure_Local[iMarker_Inlet];
      Inlet_Temperature_Total[iMarker_Inlet]                     = Inlet_Temperature_Local[iMarker_Inlet];
      Inlet_TotalTemperature_Total[iMarker_Inlet]    = Inlet_TotalTemperature_Local[iMarker_Inlet];
      Inlet_RamDrag_Total[iMarker_Inlet]                         = Inlet_RamDrag_Local[iMarker_Inlet];
      Inlet_Force_Total[iMarker_Inlet]                               = Inlet_Force_Local[iMarker_Inlet];
      Inlet_Power_Total[iMarker_Inlet]                           = Inlet_Power_Local[iMarker_Inlet];
      Inlet_Area_Total[iMarker_Inlet]                                            = Inlet_Area_Local[iMarker_Inlet];
      Inlet_XCG_Total[iMarker_Inlet]                                             = Inlet_XCG_Local[iMarker_Inlet];
      Inlet_YCG_Total[iMarker_Inlet]                                             = Inlet_YCG_Local[iMarker_Inlet];
      if (nDim == 3) Inlet_ZCG_Total[iMarker_Inlet]                                          = Inlet_ZCG_Local[iMarker_Inlet];
    }
    
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      Outlet_MassFlow_Total[iMarker_Outlet]                         = Outlet_MassFlow_Local[iMarker_Outlet];
      Outlet_Pressure_Total[iMarker_Outlet]                             = Outlet_Pressure_Local[iMarker_Outlet];
      Outlet_TotalPressure_Total[iMarker_Outlet]                = Outlet_TotalPressure_Local[iMarker_Outlet];
      Outlet_Temperature_Total[iMarker_Outlet]                  = Outlet_Temperature_Local[iMarker_Outlet];
      Outlet_TotalTemperature_Total[iMarker_Outlet] = Outlet_TotalTemperature_Local[iMarker_Outlet];
      Outlet_GrossThrust_Total[iMarker_Outlet]              = Outlet_GrossThrust_Local[iMarker_Outlet];
      Outlet_Force_Total[iMarker_Outlet]                          = Outlet_Force_Local[iMarker_Outlet];
      Outlet_Power_Total[iMarker_Outlet]                          = Outlet_Power_Local[iMarker_Outlet];
      Outlet_Area_Total[iMarker_Outlet]                                         = Outlet_Area_Local[iMarker_Outlet];
    }
    
#endif
    
    /*--- Compute the value of the average surface temperature and pressure and
     set the value in the config structure for future use ---*/
    
    for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
      if (Inlet_Area_Total[iMarker_Inlet] != 0.0) {
        Inlet_Pressure_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
        Inlet_Mach_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
        Inlet_TotalPressure_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
        Inlet_Temperature_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
        Inlet_TotalTemperature_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
        Inlet_XCG_Total[iMarker_Inlet] /= Inlet_Area_Total[iMarker_Inlet];
        Inlet_YCG_Total[iMarker_Inlet] /= Inlet_Area_Total[iMarker_Inlet];
        if (nDim == 3) Inlet_ZCG_Total[iMarker_Inlet] /= Inlet_Area_Total[iMarker_Inlet];
      }
      else {
        Inlet_Pressure_Total[iMarker_Inlet] = 0.0;
        Inlet_Mach_Total[iMarker_Inlet] = 0.0;
        Inlet_TotalPressure_Total[iMarker_Inlet] = 0.0;
        Inlet_Temperature_Total[iMarker_Inlet] = 0.0;
        Inlet_TotalTemperature_Total[iMarker_Inlet] = 0.0;
        Inlet_XCG_Total[iMarker_Inlet] = 0.0;
        Inlet_YCG_Total[iMarker_Inlet] = 0.0;
        if (nDim == 3) Inlet_ZCG_Total[iMarker_Inlet] = 0.0;
      }
      
      if (iMesh == MESH_0) {
        
        if (Engine) {
          config->SetInflow_MassFlow(iMarker_Inlet, Inlet_MassFlow_Total[iMarker_Inlet]);
          config->SetInflow_ReverseMassFlow(iMarker_Inlet, Inlet_ReverseMassFlow_Total[iMarker_Inlet]);
          config->SetInflow_Pressure(iMarker_Inlet, Inlet_Pressure_Total[iMarker_Inlet]);
          config->SetInflow_TotalPressure(iMarker_Inlet, Inlet_TotalPressure_Total[iMarker_Inlet]);
          config->SetInflow_Temperature(iMarker_Inlet, Inlet_Temperature_Total[iMarker_Inlet]);
          config->SetInflow_TotalTemperature(iMarker_Inlet, Inlet_TotalTemperature_Total[iMarker_Inlet]);
          config->SetInflow_RamDrag(iMarker_Inlet, Inlet_RamDrag_Total[iMarker_Inlet]);
          config->SetInflow_Force(iMarker_Inlet, Inlet_Force_Total[iMarker_Inlet]);
          config->SetInflow_Power(iMarker_Inlet, Inlet_Power_Total[iMarker_Inlet]);
        }
        else {
          config->SetActDiskInlet_MassFlow(iMarker_Inlet, Inlet_MassFlow_Total[iMarker_Inlet]);
          config->SetActDiskInlet_ReverseMassFlow(iMarker_Inlet, Inlet_ReverseMassFlow_Total[iMarker_Inlet]);
          config->SetActDiskInlet_Pressure(iMarker_Inlet, Inlet_Pressure_Total[iMarker_Inlet]);
          config->SetActDiskInlet_TotalPressure(iMarker_Inlet, Inlet_TotalPressure_Total[iMarker_Inlet]);
          config->SetActDiskInlet_Temperature(iMarker_Inlet, Inlet_Temperature_Total[iMarker_Inlet]);
          config->SetActDiskInlet_TotalTemperature(iMarker_Inlet, Inlet_TotalTemperature_Total[iMarker_Inlet]);
          config->SetActDiskInlet_RamDrag(iMarker_Inlet, Inlet_RamDrag_Total[iMarker_Inlet]);
          config->SetActDiskInlet_Force(iMarker_Inlet, Inlet_Force_Total[iMarker_Inlet]);
          config->SetActDiskInlet_Power(iMarker_Inlet, Inlet_Power_Total[iMarker_Inlet]);
        }
        
      }
      
    }
    
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      if (Outlet_Area_Total[iMarker_Outlet] != 0.0) {
        Outlet_Pressure_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
        Outlet_TotalPressure_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
        Outlet_Temperature_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
        Outlet_TotalTemperature_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
      }
      else {
        Outlet_Pressure_Total[iMarker_Outlet] = 0.0;
        Outlet_TotalPressure_Total[iMarker_Outlet] = 0.0;
        Outlet_Temperature_Total[iMarker_Outlet] = 0.0;
        Outlet_TotalTemperature_Total[iMarker_Outlet] = 0.0;
      }
      
      if (iMesh == MESH_0) {
        
        if (Engine) {
          config->SetExhaust_MassFlow(iMarker_Outlet, Outlet_MassFlow_Total[iMarker_Outlet]);
          config->SetExhaust_Pressure(iMarker_Outlet, Outlet_Pressure_Total[iMarker_Outlet]);
          config->SetExhaust_TotalPressure(iMarker_Outlet, Outlet_TotalPressure_Total[iMarker_Outlet]);
          config->SetExhaust_Temperature(iMarker_Outlet, Outlet_Temperature_Total[iMarker_Outlet]);
          config->SetExhaust_TotalTemperature(iMarker_Outlet, Outlet_TotalTemperature_Total[iMarker_Outlet]);
          config->SetExhaust_GrossThrust(iMarker_Outlet, Outlet_GrossThrust_Total[iMarker_Outlet]);
          config->SetExhaust_Force(iMarker_Outlet, Outlet_Force_Total[iMarker_Outlet]);
          config->SetExhaust_Power(iMarker_Outlet, Outlet_Power_Total[iMarker_Outlet]);
        }
        else {
          config->SetActDiskOutlet_MassFlow(iMarker_Outlet, Outlet_MassFlow_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_Pressure(iMarker_Outlet, Outlet_Pressure_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_TotalPressure(iMarker_Outlet, Outlet_TotalPressure_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_Temperature(iMarker_Outlet, Outlet_Temperature_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_TotalTemperature(iMarker_Outlet, Outlet_TotalTemperature_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_GrossThrust(iMarker_Outlet, Outlet_GrossThrust_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_Force(iMarker_Outlet, Outlet_Force_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_Power(iMarker_Outlet, Outlet_Power_Total[iMarker_Outlet]);
        }
        
      }
      
    }
    
    
    if (Pair) {
      
      /*--- Store delta pressure, temperature, thrust, and area ---*/
      
      for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
        
        if (Engine) {
          Inlet_TagBound = config->GetMarker_EngineInflow_TagBound(iMarker_Inlet);
          jMarker = config->GetMarker_CfgFile_EngineExhaust(Inlet_TagBound);
          Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
        }
        else {
          Inlet_TagBound = config->GetMarker_ActDiskInlet_TagBound(iMarker_Inlet);
          jMarker = config->GetMarker_CfgFile_ActDiskOutlet(Inlet_TagBound);
          Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
        }
        
        
        su2double DeltaPress = 0.0, DeltaTemp = 0.0, NetThrust = 0.0, GrossThrust = 0.0, TotalPressRatio = 0.0, TotalTempRatio = 0.0, StaticPressRatio = 0.0, StaticTempRatio = 0.0;
        
        if (Engine) {
          DeltaPress   = config->GetExhaust_Pressure(Outlet_TagBound) - config->GetInflow_Pressure(Inlet_TagBound);
          DeltaTemp         = config->GetExhaust_Temperature(Outlet_TagBound) - config->GetInflow_Temperature(Inlet_TagBound);
          NetThrust    = config->GetExhaust_GrossThrust(Outlet_TagBound) - config->GetInflow_RamDrag(Inlet_TagBound);
          GrossThrust   = config->GetExhaust_GrossThrust(Outlet_TagBound);
          TotalPressRatio   = config->GetExhaust_TotalPressure(Outlet_TagBound)/config->GetInflow_TotalPressure(Inlet_TagBound);
          TotalTempRatio    = config->GetExhaust_TotalTemperature(Outlet_TagBound)/config->GetInflow_TotalTemperature(Inlet_TagBound);
          StaticPressRatio   = config->GetExhaust_Pressure(Outlet_TagBound)/config->GetInflow_Pressure(Inlet_TagBound);
          StaticTempRatio    = config->GetExhaust_Temperature(Outlet_TagBound)/config->GetInflow_Temperature(Inlet_TagBound);
          Force = config->GetInflow_Force(Inlet_TagBound) + config->GetExhaust_Force(Outlet_TagBound);
          Power = config->GetExhaust_Power(Outlet_TagBound) - config->GetInflow_Power(Inlet_TagBound);
        }
        else {
          DeltaPress   = config->GetActDiskOutlet_Pressure(Outlet_TagBound) - config->GetActDiskInlet_Pressure(Inlet_TagBound);
          DeltaTemp         = config->GetActDiskOutlet_Temperature(Outlet_TagBound) - config->GetActDiskInlet_Temperature(Inlet_TagBound);
          NetThrust    = config->GetActDiskOutlet_GrossThrust(Outlet_TagBound) - config->GetActDiskInlet_RamDrag(Inlet_TagBound);
          GrossThrust   = config->GetActDiskOutlet_GrossThrust(Outlet_TagBound);
          TotalPressRatio   = config->GetActDiskOutlet_TotalPressure(Outlet_TagBound)/config->GetActDiskInlet_TotalPressure(Inlet_TagBound);
          TotalTempRatio    = config->GetActDiskOutlet_TotalTemperature(Outlet_TagBound)/config->GetActDiskInlet_TotalTemperature(Inlet_TagBound);
          StaticPressRatio   = config->GetActDiskOutlet_Pressure(Outlet_TagBound)/config->GetActDiskInlet_Pressure(Inlet_TagBound);
          StaticTempRatio    = config->GetActDiskOutlet_Temperature(Outlet_TagBound)/config->GetActDiskInlet_Temperature(Inlet_TagBound);
          Force = config->GetActDiskInlet_Force(Inlet_TagBound) + config->GetActDiskOutlet_Force(Outlet_TagBound);
          Power =  config->GetActDiskOutlet_Power(Outlet_TagBound) - config->GetActDiskInlet_Power(Inlet_TagBound);
          MassFlow =  config->GetActDiskInlet_MassFlow(Inlet_TagBound);
        }
        
        Mach        = Inlet_Mach_Total[iMarker_Inlet];
        Area              = Inlet_Area_Total[iMarker_Inlet];
        
        if (Engine) {
          config->SetEngine_Mach(iMarker_Inlet, Mach);
          config->SetEngine_Force(iMarker_Inlet, Force);
          config->SetEngine_Power(iMarker_Inlet, Power);
          config->SetEngine_NetThrust(iMarker_Inlet, NetThrust);
          config->SetEngine_GrossThrust(iMarker_Inlet, GrossThrust);
          config->SetEngine_Area(iMarker_Inlet, Area);
        }
        else {
          config->SetActDisk_DeltaPress(iMarker_Inlet, DeltaPress);
          config->SetActDisk_DeltaTemp(iMarker_Inlet, DeltaTemp);
          config->SetActDisk_Mach(iMarker_Inlet, Mach);
          config->SetActDisk_Force(iMarker_Inlet, Force);
          config->SetActDisk_Power(iMarker_Inlet, Power);
          config->SetActDisk_MassFlow(iMarker_Inlet, MassFlow);
          config->SetActDisk_TotalPressRatio(iMarker_Inlet, TotalPressRatio);
          config->SetActDisk_TotalTempRatio(iMarker_Inlet, TotalTempRatio);
          config->SetActDisk_StaticPressRatio(iMarker_Inlet, StaticPressRatio);
          config->SetActDisk_StaticTempRatio(iMarker_Inlet, StaticTempRatio);
          config->SetActDisk_NetThrust(iMarker_Inlet, NetThrust);
          config->SetActDisk_GrossThrust(iMarker_Inlet, GrossThrust);
          config->SetActDisk_Area(iMarker_Inlet, Area);
        }
        
      }
      
      /*--- Screen output using the values already stored in the config container ---*/
      
      if ((rank == MASTER_NODE) && (iMesh == MESH_0) ) {
        
        cout.precision(5);
        cout.setf(ios::fixed, ios::floatfield);
        
        if (write_heads && Output && !config->GetDiscrete_Adjoint()) {
          if (Engine) cout << endl   << "---------------------------- Engine properties --------------------------" << endl;
          else cout << endl   << "------------------------ Actuator Disk properties -----------------------" << endl;
        }
        
        for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
          
          if (Engine) {
            Inlet_TagBound = config->GetMarker_EngineInflow_TagBound(iMarker_Inlet);
            jMarker = config->GetMarker_CfgFile_EngineExhaust(Inlet_TagBound);
            Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
          }
          else {
            Inlet_TagBound = config->GetMarker_ActDiskInlet_TagBound(iMarker_Inlet);
            jMarker = config->GetMarker_CfgFile_ActDiskOutlet(Inlet_TagBound);
            Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
          }
          
          
          if (Engine) {
            NetThrust             =  config->GetEngine_NetThrust(iMarker_Inlet);
            GrossThrust   = config->GetEngine_GrossThrust(iMarker_Inlet);
            Power               = config->GetEngine_Power(iMarker_Inlet);
            Mach                  = config->GetEngine_Mach(iMarker_Inlet);
            Force               = config->GetEngine_Force(iMarker_Inlet);
          }
          else {
            DeltaPress      = config->GetActDisk_DeltaPress(iMarker_Inlet);
            DeltaTemp         = config->GetActDisk_DeltaTemp(iMarker_Inlet);
            TotalPressRatio       = config->GetActDisk_TotalPressRatio(iMarker_Inlet);
            TotalTempRatio        = config->GetActDisk_TotalTempRatio(iMarker_Inlet);
            StaticPressRatio      = config->GetActDisk_StaticPressRatio(iMarker_Inlet);
            StaticTempRatio           = config->GetActDisk_StaticTempRatio(iMarker_Inlet);
            NetThrust             =  config->GetActDisk_NetThrust(iMarker_Inlet);
            GrossThrust   = config->GetActDisk_GrossThrust(iMarker_Inlet);
            Power               = config->GetActDisk_Power(iMarker_Inlet);
            Mach                  = config->GetActDisk_Mach(iMarker_Inlet);
            Force               = config->GetActDisk_Force(iMarker_Inlet);
          }
          
          su2double Mach_Inf                  = config->GetMach();
          su2double Pressure_Inf  = config->GetPressure_FreeStreamND();
          
          su2double TotalPressure_Inf  = Pressure_Inf * pow( 1.0 + Mach_Inf * Mach_Inf * 0.5 * (Gamma - 1.0), Gamma     / (Gamma - 1.0));
          
          su2double MinPressure = Inlet_MinPressure_Total[iMarker_Inlet]/TotalPressure_Inf;
          su2double MaxPressure = Inlet_MaxPressure_Total[iMarker_Inlet]/TotalPressure_Inf;
          su2double AvePressure = Inlet_TotalPressure_Total[iMarker_Inlet]/TotalPressure_Inf;
          
          su2double RefDensity  = Density_Inf;
          su2double RefArea     = config->GetRefArea();
          su2double RefVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
          
          su2double Factor = (0.5*RefDensity*RefArea*RefVel2);
          su2double Ref = config->GetDensity_Ref() * config->GetVelocity_Ref() * config->GetVelocity_Ref() * 1.0 * 1.0;
          su2double DmT = GetTotal_CD() * Factor;
          
//          su2double ModDmT = 0.0;
//          if (nDim == 2) ModDmT = sqrt(GetTotal_CFx()*GetTotal_CFx() +
//                                       GetTotal_CFy()*GetTotal_CFy());
//
//          if (nDim == 3) ModDmT = sqrt(GetTotal_CFx()*GetTotal_CFx() +
//                                       GetTotal_CFy()*GetTotal_CFy() +
//                                       GetTotal_CFz()*GetTotal_CFz());
//
//          DmTVector[0] = GetTotal_CFx()/ModDmT;
//          DmTVector[1] = GetTotal_CFy()/ModDmT;
//          if (nDim == 3)  DmTVector[2] = GetTotal_CFz()/ModDmT;
          
          /*--- Set the aero drag ---*/

          su2double Aero_Drag = DmT - Force;
          su2double Aero_CD = Aero_Drag / Factor;

          SetTotal_AeroCD(Aero_CD);

          /*--- Set the solid surface drag ---*/
          
          su2double SolidSurf_Drag = DmT - Force;
          su2double SolidSurf_CD = SolidSurf_Drag / Factor;

          SetTotal_SolidCD(SolidSurf_CD);
          
          /*--- Set the net thrust value---*/
          
          su2double CT = NetThrust / Factor;

          SetTotal_NetCThrust(CT);
          
          /*--- Set the total power ---*/
          
          su2double PowerHP = Power * Ref *  config->GetVelocity_Ref() / 550.0;
          
          SetTotal_Power(PowerHP);
          
          /*--- Set the total ReverseFlow ---*/
          
          su2double ReverseFlow;
          if (Engine) ReverseFlow = fabs(config->GetInflow_ReverseMassFlow(iMarker_Inlet)  / config->GetInflow_MassFlow(Inlet_TagBound));
          else ReverseFlow = fabs(config->GetActDisk_ReverseMassFlow(iMarker_Inlet)  / config->GetActDiskInlet_MassFlow(Inlet_TagBound));
          
          SetTotal_ReverseFlow(ReverseFlow);
          
          /*--- Set the total mass flow ratio ---*/
          
          InfVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) InfVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
          if (Engine) MFR =fabs(config->GetInflow_MassFlow(Inlet_TagBound)) / (Density_Inf * sqrt(InfVel2) * config->GetHighlite_Area());
          else MFR = fabs(config->GetActDiskInlet_MassFlow(Inlet_TagBound)) / (Density_Inf * sqrt(InfVel2) * config->GetHighlite_Area());
          SetTotal_MFR(MFR);
          
          /*--- Evaluate shaft power and adiabatic efficiency (average) ---*/
          
          su2double Pstatic1, P1, P2, T1, T2;
          if (Engine) {
            Pstatic1 = config->GetInflow_Pressure(Inlet_TagBound);
            P1 = config->GetInflow_TotalPressure(Inlet_TagBound);
            P2 = config->GetExhaust_TotalPressure(Outlet_TagBound);
            T1 = config->GetInflow_TotalTemperature(Inlet_TagBound);
            T2 = config->GetExhaust_TotalTemperature(Outlet_TagBound);
          }
          else {
            Pstatic1 = config->GetActDiskInlet_Pressure(Inlet_TagBound);
            P1 = config->GetActDiskInlet_TotalPressure(Inlet_TagBound);
            P2 = config->GetActDiskOutlet_TotalPressure(Outlet_TagBound);
            T1 = config->GetActDiskInlet_TotalTemperature(Inlet_TagBound);
            T2 = config->GetActDiskOutlet_TotalTemperature(Outlet_TagBound);
          }
          
          /*-- Set the propulsive efficiency ---*/
          
          su2double mu_prop = fabs(DmT)*sqrt(RefVel2)/Power;
          SetTotal_Prop_Eff(mu_prop);
          
          /*-- Set the bypass propulsive efficiency ---*/
          
          su2double mu_bypass_prop = NetThrust*sqrt(RefVel2)/Power;
          SetTotal_ByPassProp_Eff(mu_bypass_prop);
          
          /*-- Set the fan adiabatic efficiency ---*/
          
          su2double mu_isentropic = 0.0;
          if ((P2/P1) > 0.0) mu_isentropic =    (T1/(T2-T1))*(pow((P2/P1),(Gamma-1.0)/Gamma)-1.0);
          SetTotal_Adiab_Eff(mu_isentropic);
          
          /*-- Set the polytropic efficiency ---*/
          
          su2double poly_coeff = 1.0/(1.0-log(T2/T1)/log(P2/P1));
          su2double mu_polytropic = ((Gamma-1.0)/Gamma)/((poly_coeff-1.0)/poly_coeff);
          SetTotal_Poly_Eff(mu_polytropic);
          
          if (write_heads && Output && !config->GetDiscrete_Adjoint()) {
            
            if (iMarker_Inlet > 0) cout << endl;
            
            /*--- Geometry defintion ---*/
            
            if (Engine) cout <<"Engine surfaces: " << Inlet_TagBound << ", " << Outlet_TagBound << "." << endl;
            else cout <<"Actuator disk surfaces: " << Inlet_TagBound << ", " << Outlet_TagBound << "." << endl;
            
            if (nDim == 2) {
              if (config->GetSystemMeasurements() == SI)
                cout <<"CG (m): (" << Inlet_XCG_Total[iMarker_Inlet] <<", " << Inlet_YCG_Total[iMarker_Inlet] << "). Length (m): " << Inlet_Area_Total[iMarker_Inlet] << "." << endl;
              else if (config->GetSystemMeasurements() == US)
                cout <<"CG (in): (" << Inlet_XCG_Total[iMarker_Inlet]*12.0 <<", " << Inlet_YCG_Total[iMarker_Inlet]*12.0 << "). Length (in): " << Inlet_Area_Total[iMarker_Inlet]*12.0 << "." << endl;
              cout << endl;
            }
            
            if (nDim ==3) {
              if (config->GetSystemMeasurements() == SI)
                cout <<"CG (m): (" << Inlet_XCG_Total[iMarker_Inlet] <<", " << Inlet_YCG_Total[iMarker_Inlet] <<", " << Inlet_ZCG_Total[iMarker_Inlet] << "). Area (m^2): " << Inlet_Area_Total[iMarker_Inlet] << ". Radius (m): " << sqrt(Inlet_Area_Total[iMarker_Inlet]/PI_NUMBER) << "." << endl;
              else if (config->GetSystemMeasurements() == US)
                cout <<"CG (in): (" << Inlet_XCG_Total[iMarker_Inlet]*12.0 <<", " << Inlet_YCG_Total[iMarker_Inlet]*12.0 <<", " << Inlet_ZCG_Total[iMarker_Inlet]*12.0 << "). Area (in^2): " << Inlet_Area_Total[iMarker_Inlet]*12.0*12.0 << "." << endl;
              cout << endl;
            }
            
            
            /*--- Flow field descritption  ---*/
            
            if (config->GetSystemMeasurements() == SI) {
              cout << setprecision(2) << "Inlet Ave. P (Pa): " << Pstatic1*config->GetPressure_Ref() << setprecision(3) <<  ". Inlet Ave. Mach: " << Mach << "." << endl;
              cout << setprecision(2) << "Outlet Ave. PT (Pa): " << P2*config->GetPressure_Ref() << ". Outlet Ave. TT (K): " << T2*config->GetTemperature_Ref() << "." << endl;
            }
            else if (config->GetSystemMeasurements() == US) {
              cout << setprecision(2) << "Inlet Ave. P (psf): " << Pstatic1*config->GetPressure_Ref() << setprecision(3) <<  ". Inlet Ave. Mach: " << Mach << "." << endl;
              cout << setprecision(2) << "Outlet Ave. PT (psf): " << P2*config->GetPressure_Ref() << ". Outlet Ave. TT (R): " << T2*config->GetTemperature_Ref() << "." << endl;
            }
            
            cout << "Inlet min. PT/PTinf: " << MinPressure << ". Inlet max. PT/PTinf: " << MaxPressure << ". Inlet Ave. PT/PTinf: " << AvePressure << endl;
            
            su2double InfVel2, Inlet_MassFlow, Outlet_MassFlow;
            
            if (Engine) Inlet_MassFlow = fabs(config->GetInflow_MassFlow(Inlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
            else Inlet_MassFlow = fabs(config->GetActDiskInlet_MassFlow(Inlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
            
            if (config->GetSystemMeasurements() == SI) { cout << "Inlet mass flow (kg/s): "; cout << setprecision(2) << Inlet_MassFlow; }
            else if (config->GetSystemMeasurements() == US) { cout << "Inlet mass flow (lbs/s): "; cout << setprecision(2) << Inlet_MassFlow * 32.174; }
            
            if (Engine) Outlet_MassFlow = fabs(config->GetExhaust_MassFlow(Outlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
            else Outlet_MassFlow = fabs(config->GetActDiskOutlet_MassFlow(Outlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
            
            //          if (config->GetSystemMeasurements() == SI) { cout << ". Outlet mass flow (kg/s): "; cout << setprecision(2) << Outlet_MassFlow; }
            //          else if (config->GetSystemMeasurements() == US) { cout << ". Outlet mass flow (lbs/s): "; cout << setprecision(2) << Outlet_MassFlow * 32.174; }
            
            if (Inlet_MassFlow > Outlet_MassFlow) cout << ". I/O diff.: " << setprecision(2) << 100.0*fabs(1.0-(Outlet_MassFlow/Inlet_MassFlow)) << "%";
            else cout << ". I/O diff.: " << setprecision(2) << -100.0*fabs(1.0-(Inlet_MassFlow/Outlet_MassFlow)) << "%";
            
            InfVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) InfVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
            cout << setprecision(2) << ". MFR: " << MFR << "." << endl;
            
            if (!Engine) {
              
              cout << setprecision(3) << "PT in/out ratio: " << TotalPressRatio << ". TT in/out ratio: " << TotalTempRatio  <<  "." << endl;
              
              if (config->GetActDisk_Jump() == VARIABLES_JUMP) {
                if (config->GetSystemMeasurements() == SI) cout << setprecision(3) << "P in/out jump (Pa): ";
                else if (config->GetSystemMeasurements() == US) cout << setprecision(3) << "P in/out jump (psf): ";
                cout << setprecision(3) << DeltaPress * config->GetPressure_Ref();
                if (config->GetSystemMeasurements() == SI) cout << setprecision(3) << ". T in/out jump (K): ";
                else if (config->GetSystemMeasurements() == US) cout << setprecision(3) << ". T in/out jump (R): ";
                cout << setprecision(3) << DeltaTemp * config->GetTemperature_Ref() <<"."<< endl;
              }
              else  if (config->GetActDisk_Jump() == RATIO) {
                cout << setprecision(3) << "P in/out ratio: ";
                cout << setprecision(3) << StaticPressRatio;
                cout  << setprecision(3) <<". T in/out ratio: ";
                cout << setprecision(3) << StaticTempRatio <<"."<< endl;
              }
            }
            
            cout << setprecision(1) << "\nProp. eff. (D-T.V/Shaft P): " << 100*mu_prop << "%. By-pass prop. eff. (NetT.V/Shaft P): " << 100*mu_bypass_prop <<   "%." << endl;
            cout << setprecision(1) << "Fan adiabatic eff.: " << 100*mu_isentropic  << "%. Fan poly. eff.: " << 100*mu_polytropic << "%. Poly coeff. (n): " << setprecision(4) << poly_coeff << "." << endl;
            
            cout << endl;
            
            
            /*--- Forces descritption  ---*/
            
            if (config->GetSystemMeasurements() == SI) cout << setprecision(1) << "Ram Drag (N): ";
            else if (config->GetSystemMeasurements() == US) cout << setprecision(1) << "Ram Drag (lbf): ";
            cout << (GrossThrust-NetThrust) * Ref;
            
            if (config->GetSystemMeasurements() == SI) cout << setprecision(1) << ". Gross Thrust (N): ";
            else if (config->GetSystemMeasurements() == US) cout << setprecision(1) << ". Gross Thrust (lbf): ";
            cout << -GrossThrust * Ref  << "." << endl;
            
            if (config->GetSystemMeasurements() == SI) cout << setprecision(1) << "Open surfaces Thurst (N): ";
            else if (config->GetSystemMeasurements() == US) cout  << setprecision(1) << "Open surfaces Thrust (lbf): ";
            cout<< setprecision(1) << Force * Ref << ". Open surfaces CT: " << setprecision(5) << -Force / Factor << "." << endl;
            
            if (config->GetSystemMeasurements() == SI) cout << "Solid surfaces Drag (N): ";
            else if (config->GetSystemMeasurements() == US) cout << "Solid surfaces Drag (lbf): ";
            cout << setprecision(1) << SolidSurf_Drag * Ref << ". Solid surfaces CD: " << setprecision(5) << SolidSurf_CD << "." << endl;

            if (config->GetSystemMeasurements() == SI) cout << setprecision(1) <<"Net Thrust (N): ";
            else if (config->GetSystemMeasurements() == US) cout << setprecision(1) << "Net Thrust (lbf): ";
            cout << setprecision(5) << -NetThrust * Ref  << ". Net CT: " << CT;
            
            if (config->GetSystemMeasurements() == SI) {
              cout << ". Power (W): ";
              cout << setprecision(1) << Power * Ref *  config->GetVelocity_Ref()  << "." << endl;
            }
            else if (config->GetSystemMeasurements() == US) {
              cout << ". Power (HP): ";
              cout << setprecision(1) << Power * Ref *  config->GetVelocity_Ref() / 550.0 << "." << endl;
            }
            
          }
          
        }
        
        if (write_heads && Output && !config->GetDiscrete_Adjoint()) cout << "-------------------------------------------------------------------------" << endl << endl;
        
      }
      
    }
    
    
    delete [] Outlet_MassFlow_Local;
    delete [] Outlet_Temperature_Local;
    delete [] Outlet_TotalTemperature_Local;
    delete [] Outlet_Pressure_Local;
    delete [] Outlet_TotalPressure_Local;
    delete [] Outlet_Area_Local;
    delete [] Outlet_GrossThrust_Local;
    delete [] Outlet_Force_Local;
    delete [] Outlet_Power_Local;
    
    delete [] Outlet_MassFlow_Total;
    delete [] Outlet_Temperature_Total;
    delete [] Outlet_TotalTemperature_Total;
    delete [] Outlet_Pressure_Total;
    delete [] Outlet_TotalPressure_Total;
    delete [] Outlet_Area_Total;
    delete [] Outlet_GrossThrust_Total;
    delete [] Outlet_Force_Total;
    delete [] Outlet_Power_Total;
    
    delete [] Inlet_MassFlow_Local;
    delete [] Inlet_ReverseMassFlow_Local;
    delete [] Inlet_Temperature_Local;
    delete [] Inlet_TotalTemperature_Local;
    delete [] Inlet_Pressure_Local;
    delete [] Inlet_Mach_Local;
    delete [] Inlet_MinPressure_Local;
    delete [] Inlet_MaxPressure_Local;
    delete [] Inlet_TotalPressure_Local;
    delete [] Inlet_Area_Local;
    delete [] Inlet_RamDrag_Local;
    delete [] Inlet_Force_Local;
    delete [] Inlet_Power_Local;
    delete [] Inlet_XCG_Local;
    delete [] Inlet_YCG_Local;
    delete [] Inlet_ZCG_Local;
    
    delete [] Inlet_MassFlow_Total;
    delete [] Inlet_ReverseMassFlow_Total;
    delete [] Inlet_Temperature_Total;
    delete [] Inlet_TotalTemperature_Total;
    delete [] Inlet_Pressure_Total;
    delete [] Inlet_Mach_Total;
    delete [] Inlet_MinPressure_Total;
    delete [] Inlet_MaxPressure_Total;
    delete [] Inlet_TotalPressure_Total;
    delete [] Inlet_Area_Total;
    delete [] Inlet_RamDrag_Total;
    delete [] Inlet_Force_Total;
    delete [] Inlet_Power_Total;
    delete [] Inlet_XCG_Total;
    delete [] Inlet_YCG_Total;
    delete [] Inlet_ZCG_Total;
    
    delete [] Inlet_MassFlow;
    delete [] Inlet_Mach;
    delete [] Inlet_MinPressure;
    delete [] Inlet_MaxPressure;
    delete [] Inlet_ReverseMassFlow;
    delete [] Inlet_Pressure;
    delete [] Inlet_TotalPressure;
    delete [] Inlet_Temperature;
    delete [] Inlet_TotalTemperature;
    delete [] Inlet_Area;
    delete [] Inlet_RamDrag;
    delete [] Inlet_Force;
    delete [] Inlet_Power;
    delete [] Inlet_XCG;
    delete [] Inlet_YCG;
    delete [] Inlet_ZCG;
    
    delete [] Outlet_MassFlow;
    delete [] Outlet_Pressure;
    delete [] Outlet_TotalPressure;
    delete [] Outlet_Temperature;
    delete [] Outlet_TotalTemperature;
    delete [] Outlet_Area;
    delete [] Outlet_GrossThrust;
    delete [] Outlet_Force;
    delete [] Outlet_Power;
    
  }
  
}

void CEulerSolver::SetActDisk_BCThrust(CGeometry *geometry, CSolver **solver_container,
                                       CConfig *config, unsigned short iMesh, bool Output) {
  
  su2double Massflow = 0.0 , Target_Massflow = 0.0, DragMinusThrust = 0.0 , Target_DragMinusThrust = 0.0, Target_NetThrust = 0.0, BCThrust = 0.0, BCThrust_inc = 0.0;
  unsigned short iDim, iMarker;
  unsigned long iVertex, iPoint;
  su2double  *V_inlet = NULL, Pressure,
  Density, T0_Ti, ATerm, BTerm, LHS, RHS, RHS_PDelta, RHS_MDelta, F, DF_DLa, CTerm_, DTerm_,
  ETerm, La, La_old, TotalArea, To_Ti, DeltaT, Po_Pi, DeltaP, Area, Velocity_Normal,
  SoundSpeed2, Force_Normal,
  RefDensity, RefArea, RefVel2, Factor, Ref;
  unsigned short iter;
  string Marker_Tag;
  su2double Target_Force, Force, Target_Power, Power, NetThrust, BCThrust_old, Initial_BCThrust;
  bool ActDisk_Info;
  su2double MyBCThrust, BCThrust_Init;
  
  su2double dNetThrust_dBCThrust        = config->GetdNetThrust_dBCThrust();
  unsigned short Kind_ActDisk           = config->GetKind_ActDisk();
  bool ratio                            = (config->GetActDisk_Jump() == RATIO);
  unsigned long Update_BCThrust         = config->GetUpdate_BCThrust();
  unsigned long Iter_Fixed_NetThrust    = config->GetIter_Fixed_NetThrust();
  unsigned long ExtIter                 = config->GetExtIter();
  bool Update_BCThrust_Bool             = false;
  bool restart                          = (config->GetRestart() || config->GetRestart_Flow());
  su2double Fan_Poly_Eff                = config->GetFan_Poly_Eff();
  su2double PolyCoeff                   = 1.0/(1.0-((Gamma-1.0)/Gamma)/Fan_Poly_Eff);
  
  RefDensity   = Density_Inf;
  RefArea = config->GetRefArea();
  RefVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  
  Factor = (0.5*RefDensity*RefArea*RefVel2);
  Ref = config->GetDensity_Ref() * config->GetVelocity_Ref() * config->GetVelocity_Ref() * 1.0 * 1.0;
  
  /*--- Delta P and delta T are inputs ---*/
  
  if (Kind_ActDisk == VARIABLES_JUMP) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
        
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        
        if (ratio) {
          if (config->GetMach()  < 0.5) {
            DeltaP       = config->GetActDisk_PressJump(Marker_Tag, 0);
            DeltaT       = config->GetActDisk_TempJump(Marker_Tag, 0);
          }
          else {
            DeltaP       = config->GetActDisk_PressJump(Marker_Tag, 1);
            DeltaT       = config->GetActDisk_TempJump(Marker_Tag, 1);
          }
        }
        else {
          if (config->GetMach()  < 0.5) {
            DeltaP       = max(0.0, config->GetActDisk_PressJump(Marker_Tag, 0) / config->GetPressure_Ref());
            DeltaT       = max(0.0, config->GetActDisk_TempJump(Marker_Tag, 0) / config->GetTemperature_Ref());
          }
          else {
            DeltaP       = max(0.0, config->GetActDisk_PressJump(Marker_Tag, 1) / config->GetPressure_Ref());
            DeltaT       = max(0.0, config->GetActDisk_TempJump(Marker_Tag, 1) / config->GetTemperature_Ref());
          }
        }
        
        /*--- Set the Delta P, Delta T values at each discrete point (uniform distribution)  ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          SetActDisk_DeltaP(iMarker, iVertex, DeltaP);
          SetActDisk_DeltaT(iMarker, iVertex, DeltaT);
        }
        
      }
    }
  }
  
  /*--- Iteration using BCThrust ---*/
  
  else {
    
    if (ExtIter == 0) BCThrust_Counter = 0;
    
    /*--- Only the fine mesh level should check the convergence criteria ---*/
    
    if ((iMesh == MESH_0) && Output) {
      
      /*--- Initialize the update flag to false ---*/
      
      Update_BCThrust_Bool = false;
      
      /*--- Reevaluate BCThrust at a fix number of iterations ---*/
      
      if ((ExtIter % Iter_Fixed_NetThrust == 0) && (ExtIter != 0)) {
        BCThrust_Counter++;
        if ((BCThrust_Counter != 0) &&
            (BCThrust_Counter != 1) &&
            (BCThrust_Counter != Update_BCThrust) &&
            (BCThrust_Counter != Update_BCThrust + 2) &&
            (BCThrust_Counter != Update_BCThrust + 4) ) Update_BCThrust_Bool = true;
        else Update_BCThrust_Bool = false;
      }
      
      /*--- Store the update boolean for use on other mesh levels in the MG ---*/
      
      config->SetUpdate_BCThrust_Bool(Update_BCThrust_Bool);
      
    }
    
    else {
      Update_BCThrust_Bool = config->GetUpdate_BCThrust_Bool();
    }
    
    
    /*--- If it is the first iteration, set the BCThrust to a meaning full target value,
     * this can be done at an initialization level, for the time being it is OK here ---*/
    
    if (ExtIter == 0) {
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
            (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          
          if (Kind_ActDisk == NET_THRUST) {
            if (restart)
              Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
            else {
              if (config->GetMach() < 0.5) Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
              else Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
            }
            config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
            config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
          }
          
          if (Kind_ActDisk == BC_THRUST) {
            if (restart)
              Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
            else {
              if (config->GetMach() < 0.5) Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
              else Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
            }
            config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
            config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
          }
          
          if (Kind_ActDisk == POWER) {
            Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
            config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
            config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
          }
          
          if (Kind_ActDisk == DRAG_MINUS_THRUST) {
            Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
            config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
            config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
          }
          
          if (Kind_ActDisk == MASSFLOW) {
            Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
            config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
            config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
          }
          
        }
      }
    }
    
    /*--- Typical iteration to set the value of BC Thrust at each actuator disk ---*/
    
    if (Update_BCThrust_Bool && Output) {
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        
        if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
            (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
          
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          
          if (Kind_ActDisk == NET_THRUST) {
            
            if (config->GetMach() < 0.5) {
              Target_NetThrust    = fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
            }
            else {
              Target_NetThrust    = fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
            }
            NetThrust    = config->GetActDisk_NetThrust(Marker_Tag);
            BCThrust_old = config->GetActDisk_BCThrust_Old(Marker_Tag);
            BCThrust_inc = (1.0/dNetThrust_dBCThrust)*(Target_NetThrust - NetThrust);

            if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
            else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            if (iMesh == MESH_0) {
              config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
              BCThrust_Init = BCThrust*Ref;
              config->SetInitial_BCThrust(BCThrust_Init);
            }
          }
          
          if (Kind_ActDisk == BC_THRUST) {
            
            if (config->GetMach() < 0.5) {
              Target_Force =  fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
            }
            else {
              Target_Force =  fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
            }
            
            Force        = -config->GetActDisk_Force(Marker_Tag);
            BCThrust_old = config->GetActDisk_BCThrust_Old(Marker_Tag);
            BCThrust_inc = (1.0/dNetThrust_dBCThrust)*(Target_Force - Force);
            if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
            else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            
            if (iMesh == MESH_0) {
              config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
              BCThrust_Init = BCThrust*Ref;
              config->SetInitial_BCThrust(BCThrust_Init);
            }
            
          }
          
          if (Kind_ActDisk == POWER) {
            
            if (config->GetMach() < 0.5) {
              Target_Power =  fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / (Ref * config->GetVelocity_Ref() /  550.0));
            }
            else {
              Target_Power =  fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / (Ref * config->GetVelocity_Ref() /  550.0));
            }
            
            Power        = config->GetActDisk_Power(Marker_Tag);
            BCThrust_old = config->GetActDisk_BCThrust_Old(Marker_Tag);
            BCThrust_inc = (1.0/dNetThrust_dBCThrust)*(Target_Power - Power);
            if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
            else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            if (iMesh == MESH_0) {
              config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
              BCThrust_Init = BCThrust*Ref;
              config->SetInitial_BCThrust(BCThrust_Init);
            }
          }
          
          if (Kind_ActDisk == DRAG_MINUS_THRUST) {
            
            if (config->GetMach() < 0.5) {
              Target_DragMinusThrust  =  -fabs(config->GetActDisk_PressJump(Marker_Tag, 0)) * Factor;
            }
            else {
              Target_DragMinusThrust  =  -fabs(config->GetActDisk_PressJump(Marker_Tag, 1)) * Factor;
            }
            
            DragMinusThrust = GetTotal_CD() * Factor;
            BCThrust_old    = config->GetActDisk_BCThrust_Old(Marker_Tag);
            BCThrust_inc    = -(1.0/dNetThrust_dBCThrust)*(Target_DragMinusThrust - DragMinusThrust);
            if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
            else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            if (iMesh == MESH_0) {
              config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
              BCThrust_Init = BCThrust*Ref;
              config->SetInitial_BCThrust(BCThrust_Init);
            }
            
          }
          
          if (Kind_ActDisk == MASSFLOW) {
            
            if (config->GetMach() < 0.5) {
              Target_Massflow  =  fabs(config->GetActDisk_PressJump(Marker_Tag, 0) / (config->GetDensity_Ref() * config->GetVelocity_Ref()));
              if (config->GetSystemMeasurements() == US) Target_Massflow /= 32.174;
            }
            else {
              Target_Massflow  =  fabs(config->GetActDisk_PressJump(Marker_Tag, 1) / (config->GetDensity_Ref() * config->GetVelocity_Ref()));
              if (config->GetSystemMeasurements() == US) Target_Massflow /= 32.174;
            }
            
            Massflow = config->GetActDisk_MassFlow(Marker_Tag);
            BCThrust_old    = config->GetActDisk_BCThrust_Old(Marker_Tag);
            BCThrust_inc    = (1.0/dNetThrust_dBCThrust)*(Target_Massflow - Massflow);
            if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
            else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            if (iMesh == MESH_0) {
              config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
              BCThrust_Init = BCThrust*Ref;
              config->SetInitial_BCThrust(BCThrust_Init);
            }
            
          }
          
        }
        
      }
      
      /*--- After a complete update of BC_Thrust
       update the value of BC Thrust (old) for future iterations ---*/
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
            (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if ((Kind_ActDisk == NET_THRUST) || (Kind_ActDisk == BC_THRUST) ||
              (Kind_ActDisk == POWER) || (Kind_ActDisk == DRAG_MINUS_THRUST) ||
              (Kind_ActDisk == MASSFLOW)) {
            BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            config->SetActDisk_BCThrust_Old(Marker_Tag, BCThrust);
          }
        }
      }
      
    }
    
    /*--- Evaluate the pressure jump at each node using the total thrust ---*/
    
    if ((Update_BCThrust_Bool && Output) || (ExtIter == 0)) {
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        
        if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
            (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
          
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          RefDensity  = Density_Inf;
          RefArea = config->GetRefArea();
          RefVel2 = 0.0; for (iDim = 0; iDim < nDim; iDim++) RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
          
          Factor = (0.5*RefDensity*RefArea*RefVel2);
          Ref = config->GetDensity_Ref() * config->GetVelocity_Ref() * config->GetVelocity_Ref() * 1.0 * 1.0;
          BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            
            if (geometry->node[iPoint]->GetDomain()) {
              
              geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
              
              if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) {
                for (iDim = 0; iDim < nDim; iDim++) { Vector[iDim] = -Vector[iDim]; }
              }
              
              Area = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) { Area += Vector[iDim]*Vector[iDim]; }
              Area = sqrt (Area);
              
              /*--- Use the inlet state to compute the Pressure and Temperature jumps ---*/
              
              if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET)
                V_inlet = node[iPoint]->GetPrimitive();
              if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)
                V_inlet = GetDonorPrimVar(iMarker, iVertex);
              
              Density       = V_inlet[nDim+2];
              Pressure      = V_inlet[nDim+1];
              SoundSpeed2   = Pressure*Gamma/Density;
              TotalArea     = config->GetActDisk_Area(Marker_Tag);
              Force_Normal  = Area*(BCThrust/TotalArea);
              
              
              Velocity_Normal = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) {
                Velocity_Normal += V_inlet[iDim+1]*Vector[iDim]/Area;
              }
              
              
              if (Velocity_Normal > EPS) {
                
                /*--- Ratio of the total temperature to the temperature at the inflow ---*/
                
                T0_Ti = 1.0 + ((Gamma-1.0)/SoundSpeed2)*(0.5*Velocity_Normal*Velocity_Normal + Force_Normal/(Density*Area));
                
                
                ATerm = 2.0*T0_Ti/(Gamma+1.0);
                BTerm = 0.5*(Gamma+1.0)/(Gamma-1.0);
                LHS = fabs(Velocity_Normal)/(sqrt(SoundSpeed2)*pow(ATerm,BTerm));
                
                CTerm_ = (PolyCoeff-1.0)/(PolyCoeff+1.0);
                DTerm_ = 1.0/(PolyCoeff-1.0);
                
                La = EPS; La_old = EPS;
                
                for (iter = 0; iter < 100; iter++) {
                  
                  ETerm = ((1.0-CTerm_*La*La)/(1.0-CTerm_+EPS));
                  
                  RHS = La*pow(ETerm, DTerm_);
                  
                  ETerm = ((1.0-CTerm_*(La+1E-6)*(La+1E-6))/(1.0-CTerm_+EPS));
                  RHS_PDelta = (La+1E-6)*pow(ETerm, DTerm_);
                  
                  ETerm = ((1.0-CTerm_*(La-1E-6)*(La-1E-6))/(1.0-CTerm_+EPS));
                  RHS_MDelta = (La-1E-6)*pow(ETerm, DTerm_);
                  
                  /*--- Objective function and finitte differences derivative ---*/
                  
                  F = RHS - LHS;
                  DF_DLa = (RHS_PDelta - RHS_MDelta)/2E-6;
                  
                  /*--- Newton's step ---*/
                  
                  La_old = La;
                  La = La_old - 0.75*(F/DF_DLa);
                  
                  if (fabs(F) < 1E-10) break;
                  
                }
                
                if (iter == 99) cout << "The laval number evaluation is not converging." << endl;
                
                /*--- Laval is bounded ---*/
                
                La = min(La, sqrt(6.0));  La = max(La, 0.0);
                
                To_Ti = max(1.0, T0_Ti*(1.0-CTerm_*La*La));
                SetActDisk_DeltaT(iMarker, iVertex, To_Ti);
                
                Po_Pi = max(1.0, pow(To_Ti, PolyCoeff*DTerm_));
                SetActDisk_DeltaP(iMarker, iVertex, Po_Pi);
                
              }
              else {
                SetActDisk_DeltaT(iMarker, iVertex, 1.0);
                SetActDisk_DeltaP(iMarker, iVertex, 1.0);
              }
              
            }
            
          }
        }
        
      }
    }
  }
  
  /*--- Broadcast some information to the master node ---*/
  
  ActDisk_Info = false;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
        (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
      ActDisk_Info = true;
    }
  }
  if (!ActDisk_Info) config->SetInitial_BCThrust(0.0);
  
  MyBCThrust = config->GetInitial_BCThrust();
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&MyBCThrust, &BCThrust, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  BCThrust = MyBCThrust;
#endif
  config->SetInitial_BCThrust(BCThrust);
  
}

void CEulerSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, bool Output) {
  
  su2double Target_CL = 0.0, AoA = 0.0, Vel_Infty[3], AoA_inc = 0.0, Vel_Infty_Mag, Old_AoA,
  dCL_dAlpha_, dCD_dCL_, dCMx_dCL_, dCMy_dCL_, dCMz_dCL_;
  unsigned long Wrt_Con_Freq;
  unsigned short iDim;
  
  unsigned long Iter_Fixed_CL = config->GetIter_Fixed_CL();
  unsigned long Update_Alpha = config->GetUpdate_Alpha();
  
  unsigned long ExtIter       = config->GetExtIter();
  bool write_heads = ((ExtIter % Iter_Fixed_CL == 0) && (ExtIter != 0));
  su2double Beta                 = config->GetAoS()*PI_NUMBER/180.0;
  su2double dCL_dAlpha           = config->GetdCL_dAlpha()*180.0/PI_NUMBER;
  bool Update_AoA             = false;
  
  if (ExtIter == 0) AoA_Counter = 0;
  
  /*--- Only the fine mesh level should check the convergence criteria ---*/
  
  if ((iMesh == MESH_0) && Output) {
    
    /*--- Initialize the update flag to false ---*/
    
    Update_AoA = false;

    /*--- Reevaluate Angle of Attack at a fixed number of iterations ---*/
    
    if ((ExtIter % Iter_Fixed_CL == 0) && (ExtIter != 0)) {
      AoA_Counter++;
      if ((AoA_Counter <= Update_Alpha)) Update_AoA = true;
      else Update_AoA = false;
    }
    
    /*--- Store the update boolean for use on other mesh levels in the MG ---*/
    
    config->SetUpdate_AoA(Update_AoA);
    
  }
  
  else {
    Update_AoA = config->GetUpdate_AoA();
  }
  
  if (Update_AoA && Output) {
    
    /*--- Retrieve the specified target CL value. ---*/
    
    Target_CL = config->GetTarget_CL();
    
    /*--- Retrieve the old AoA (radians) ---*/
    
    AoA_old = config->GetAoA()*PI_NUMBER/180.0;
    
    /*--- Estimate the increment in AoA based on dCL_dAlpha (radians) ---*/
    
    AoA_inc = (1.0/dCL_dAlpha)*(Target_CL - Total_CL);
    
    /*--- Compute a new value for AoA on the fine mesh only (radians)---*/
    
    if (iMesh == MESH_0) AoA = AoA_old + AoA_inc;
    else { AoA = config->GetAoA()*PI_NUMBER/180.0; }
    
    /*--- Only the fine mesh stores the updated values for AoA in config ---*/
    
    if (iMesh == MESH_0) {
      config->SetAoA(AoA*180.0/PI_NUMBER);
    }
    
    /*--- Update the freestream velocity vector at the farfield ---*/
    
    for (iDim = 0; iDim < nDim; iDim++)
      Vel_Infty[iDim] = GetVelocity_Inf(iDim);
    
    /*--- Compute the magnitude of the free stream velocity ---*/
    
    Vel_Infty_Mag = 0;
    for (iDim = 0; iDim < nDim; iDim++)
      Vel_Infty_Mag += Vel_Infty[iDim]*Vel_Infty[iDim];
    Vel_Infty_Mag = sqrt(Vel_Infty_Mag);
    
    /*--- Compute the new freestream velocity with the updated AoA ---*/
    
    if (nDim == 2) {
      Vel_Infty[0] = cos(AoA)*Vel_Infty_Mag;
      Vel_Infty[1] = sin(AoA)*Vel_Infty_Mag;
    }
    if (nDim == 3) {
      Vel_Infty[0] = cos(AoA)*cos(Beta)*Vel_Infty_Mag;
      Vel_Infty[1] = sin(Beta)*Vel_Infty_Mag;
      Vel_Infty[2] = sin(AoA)*cos(Beta)*Vel_Infty_Mag;
    }
    
    /*--- Store the new freestream velocity vector for the next iteration ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_Inf[iDim] = Vel_Infty[iDim];
    }
    
    /*--- Only the fine mesh stores the updated values for velocity in config ---*/
    
    if (iMesh == MESH_0) {
      for (iDim = 0; iDim < nDim; iDim++)
        config->SetVelocity_FreeStreamND(Vel_Infty[iDim], iDim);
    }
    
    /*--- Output some information to the console with the headers ---*/
    
    if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && !config->GetDiscrete_Adjoint()) {
      Old_AoA = config->GetAoA() - AoA_inc*(180.0/PI_NUMBER);
      
      cout.precision(7);
      cout.setf(ios::fixed, ios::floatfield);
      cout << endl << "----------------------------- Fixed CL Mode -----------------------------" << endl;
      cout << "CL: " << Total_CL;
      cout << " (target: " << config->GetTarget_CL() <<")." << endl;
      cout.precision(4);
      cout << "Previous AoA: " << Old_AoA << " deg";
      cout << ", new AoA: " << config->GetAoA() << " deg." << endl;
      cout << "-------------------------------------------------------------------------" << endl << endl;
    }

  }

	unsigned long Iter_dCL_dAlpha = config->GetIter_dCL_dAlpha();

  if ((config->GetnExtIter()-Iter_dCL_dAlpha == ExtIter) && Output) {

    AoA_old = config->GetAoA();

    if (config->GetnExtIter()-Iter_dCL_dAlpha == ExtIter) {
      Wrt_Con_Freq = SU2_TYPE::Int(su2double(config->GetIter_dCL_dAlpha())/10.0);
      config->SetWrt_Con_Freq(Wrt_Con_Freq);
      Total_CD_Prev = Total_CD;
      Total_CL_Prev = Total_CL;
      Total_CMx_Prev = Total_CMx;
      Total_CMy_Prev = Total_CMy;
      Total_CMz_Prev = Total_CMz;
      AoA_inc = 0.001;
      AoA_FD_Change = true;
    }
    
    if ((rank == MASTER_NODE) && (iMesh == MESH_0) && !config->GetDiscrete_Adjoint()) {
      
    	if (config->GetnExtIter()-Iter_dCL_dAlpha == ExtIter) {
       cout << endl << "----------------------------- Fixed CL Mode -----------------------------" << endl;
       cout << " Change AoA by +0.001 deg to evaluate gradient." << endl;
       cout << "-------------------------------------------------------------------------" << endl << endl;
    	}

    }

    /*--- Compute a new value for AoA on the fine mesh only (radians)---*/

    if (iMesh == MESH_0) AoA = AoA_old + AoA_inc;
    else { AoA = config->GetAoA(); }

    /*--- Only the fine mesh stores the updated values for AoA in config ---*/

    if (iMesh == MESH_0) { config->SetAoA(AoA); }

    /*--- Update the freestream velocity vector at the farfield ---*/

    for (iDim = 0; iDim < nDim; iDim++)
      Vel_Infty[iDim] = GetVelocity_Inf(iDim);

    /*--- Compute the magnitude of the free stream velocity ---*/

    Vel_Infty_Mag = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Vel_Infty_Mag += Vel_Infty[iDim]*Vel_Infty[iDim];
    Vel_Infty_Mag = sqrt(Vel_Infty_Mag);

    /*--- Compute the new freestream velocity with the updated AoA ---*/

    if (nDim == 2) {
      Vel_Infty[0] = cos(AoA*PI_NUMBER/180.0)*Vel_Infty_Mag;
      Vel_Infty[1] = sin(AoA*PI_NUMBER/180.0)*Vel_Infty_Mag;
    }
    if (nDim == 3) {
      Vel_Infty[0] = cos(AoA*PI_NUMBER/180.0)*cos(Beta)*Vel_Infty_Mag;
      Vel_Infty[1] = sin(Beta)*Vel_Infty_Mag;
      Vel_Infty[2] = sin(AoA*PI_NUMBER/180.0)*cos(Beta)*Vel_Infty_Mag;
    }

    /*--- Store the new freestream velocity vector for the next iteration ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_Inf[iDim] = Vel_Infty[iDim];
    }

    /*--- Only the fine mesh stores the updated values for velocity in config ---*/

    if (iMesh == MESH_0) {
      for (iDim = 0; iDim < nDim; iDim++)
        config->SetVelocity_FreeStreamND(Vel_Infty[iDim], iDim);
    }

  }
  
  if (AoA_FD_Change && (config->GetnExtIter()-1 == ExtIter) && Output && (iMesh == MESH_0) && !config->GetDiscrete_Adjoint()) {

    /*--- Update angle of attack ---*/

    AoA_old = config->GetAoA();
    AoA = AoA_old - 0.001;
    config->SetAoA(AoA);
    
    /*--- Use finite differences to compute  ---*/

  		dCL_dAlpha_ = (Total_CL-Total_CL_Prev)/0.001;
  		dCD_dCL_    = (Total_CD-Total_CD_Prev)/(Total_CL-Total_CL_Prev);
  		dCMx_dCL_   = (Total_CMx-Total_CMx_Prev)/(Total_CL-Total_CL_Prev);
  		dCMy_dCL_   = (Total_CMy-Total_CMy_Prev)/(Total_CL-Total_CL_Prev);
  		dCMz_dCL_   = (Total_CMz-Total_CMz_Prev)/(Total_CL-Total_CL_Prev);

  		/*--- Set the value of the  dOF/dCL in the config file ---*/

  		config->SetdCD_dCL(dCD_dCL_);
  		config->SetdCMx_dCL(dCMx_dCL_);
  		config->SetdCMy_dCL(dCMy_dCL_);
  		config->SetdCMz_dCL(dCMz_dCL_);

  		config->SetdCL_dAlpha(dCL_dAlpha_);

    if (rank == MASTER_NODE) {
      cout << endl << "----------------------------- Fixed CL Mode -----------------------------" << endl;
      cout << "Approx. Delta CL / Delta AoA: " << dCL_dAlpha_ << " (1/deg)." << endl;
      cout << "Approx. Delta CD / Delta CL: " << dCD_dCL_ << ". " << endl;
      if (nDim == 3 ) {
        cout << "Approx. Delta CMx / Delta CL: " << dCMx_dCL_ << ". " << endl;
        cout << "Approx. Delta CMy / Delta CL: " << dCMy_dCL_ << ". " << endl;
      }
      cout << "Approx. Delta CMz / Delta CL: " << dCMz_dCL_ << ". " << endl;
      cout << "-------------------------------------------------------------------------" << endl << endl;
    }
    
  }

}

void CEulerSolver::SetInlet(CConfig *config) {

  unsigned short iMarker, iDim;
  unsigned long iVertex;
  string Marker_Tag;

  /* --- Initialize quantities for inlet boundary
   * This routine does not check if the inlet boundaries are set to custom
   * values (available through the py wrapper). This is intentional; the
   * default values for these custom BCs are initialized with the default
   * values specified in the config (avoiding non physical values) --- */
  for(iMarker=0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
      string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
      su2double p_total = config->GetInlet_Ptotal(Marker_Tag);
      su2double t_total = config->GetInlet_Ttotal(Marker_Tag);
      su2double* flow_dir = config->GetInlet_FlowDir(Marker_Tag);

      for(iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
        Inlet_Ttotal[iMarker][iVertex] = t_total;
        Inlet_Ptotal[iMarker][iVertex] = p_total;
        for (iDim = 0; iDim < nDim; iDim++)
          Inlet_FlowDir[iMarker][iVertex][iDim] = flow_dir[iDim];
      }
    }
  }

}

void CEulerSolver::UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config){

  unsigned short nMGlevel, iMarker;

  // TODO: Update the fluid boundary conditions for MG
  nMGlevel = config->GetnMGLevels();
  if (nMGlevel > 1) {
    for (iMarker=0; iMarker < nMarker; iMarker++) {
      bool isCustomizable = config->GetMarker_All_PyCustom(iMarker);
      bool isInlet = (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW);
      if (isCustomizable && isInlet)
        SU2_MPI::Error("Custom inlet BCs are not currently compatible with multigrid.", CURRENT_FUNCTION);
    }
  }
}

void CEulerSolver::Evaluate_ObjFunc(CConfig *config) {
  
  unsigned short iMarker_Monitoring, Kind_ObjFunc;
  su2double Weight_ObjFunc;
  
  Total_ComboObj = 0.0;

  /*--- Loop over all monitored markers, add to the 'combo' objective ---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {

    Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);
    Kind_ObjFunc = config->GetKind_ObjFunc(iMarker_Monitoring);

    switch(Kind_ObjFunc) {
      case DRAG_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CD[iMarker_Monitoring]);
        if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCD_dCL()*(Surface_CL[iMarker_Monitoring]);
        if (config->GetFixed_CM_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCD_dCMy()*(Surface_CMy[iMarker_Monitoring]);
        break;
      case LIFT_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CL[iMarker_Monitoring]);
        break;
      case SIDEFORCE_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CSF[iMarker_Monitoring]);
        break;
      case EFFICIENCY:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CEff[iMarker_Monitoring]);
        break;
      case MOMENT_X_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CMx[iMarker_Monitoring]);
        if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCMx_dCL()*(Surface_CL[iMarker_Monitoring]);
        break;
      case MOMENT_Y_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CMy[iMarker_Monitoring]);
        if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCMy_dCL()*(Surface_CL[iMarker_Monitoring]);
        break;
      case MOMENT_Z_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CMz[iMarker_Monitoring]);
        if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCMz_dCL()*(Surface_CL[iMarker_Monitoring]);
        break;
      case FORCE_X_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*Surface_CFx[iMarker_Monitoring];
        break;
      case FORCE_Y_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*Surface_CFy[iMarker_Monitoring];
        break;
      case FORCE_Z_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*Surface_CFz[iMarker_Monitoring];
        break;
      case TOTAL_HEATFLUX:
        Total_ComboObj+=Weight_ObjFunc*Surface_HF_Visc[iMarker_Monitoring];
        break;
      case MAXIMUM_HEATFLUX:
        Total_ComboObj+=Weight_ObjFunc*Surface_MaxHF_Visc[iMarker_Monitoring];
        break;
        
        /*--- The following are not per-surface, and as a result will be
         double-counted iff multiple surfaces are specified as well as multi-objective ---*/
        
      case EQUIVALENT_AREA:
        Total_ComboObj+=Weight_ObjFunc*Total_CEquivArea;
        break;
      case NEARFIELD_PRESSURE:
        Total_ComboObj+=Weight_ObjFunc*Total_CNearFieldOF;
        break;
      case INVERSE_DESIGN_PRESSURE:
        Total_ComboObj+=Weight_ObjFunc*Total_CpDiff;
        break;
      case INVERSE_DESIGN_HEATFLUX:
        Total_ComboObj+=Weight_ObjFunc*Total_HeatFluxDiff;
        break;
      case THRUST_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*Total_CT;
        break;
      case TORQUE_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*Total_CQ;
        break;
      case FIGURE_OF_MERIT:
        Total_ComboObj+=Weight_ObjFunc*Total_CMerit;
        break;
      case SURFACE_TOTAL_PRESSURE:
        Total_ComboObj+=Weight_ObjFunc*config->GetSurface_TotalPressure(0);
        break;
      case SURFACE_STATIC_PRESSURE:
        Total_ComboObj+=Weight_ObjFunc*config->GetSurface_Pressure(0);
        break;
      case SURFACE_MASSFLOW:
        Total_ComboObj+=Weight_ObjFunc*config->GetSurface_MassFlow(0);
        break;
      case SURFACE_MACH:
        Total_ComboObj+=Weight_ObjFunc*config->GetSurface_Mach(0);
        break;
      case CUSTOM_OBJFUNC:
        Total_ComboObj+=Weight_ObjFunc*Total_Custom_ObjFunc;
        break;
      default:
        break;

    }
  }
  
}

void CEulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, iVar, jVar, kVar, jDim;
  unsigned long iPoint, iVertex;
  su2double *Normal = NULL, *GridVel = NULL, Area, UnitNormal[3], *NormalArea,
  ProjGridVel = 0.0, turb_ke;
  su2double Density_b, StaticEnergy_b, Enthalpy_b, *Velocity_b, Kappa_b, Chi_b, Energy_b, VelMagnitude2_b, Pressure_b;
  su2double Density_i, *Velocity_i, ProjVelocity_i = 0.0, Energy_i, VelMagnitude2_i;
  su2double **Jacobian_b, **DubDu;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  
  Normal = new su2double[nDim];
  NormalArea = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  Velocity_i = new su2double[nDim];
  Jacobian_b = new su2double*[nVar];
  DubDu = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_b[iVar] = new su2double[nVar];
    DubDu[iVar] = new su2double[nVar];
  }
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negative for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++) {
        NormalArea[iDim] = -Normal[iDim];
        UnitNormal[iDim] = -Normal[iDim]/Area;
      }
      
      /*--- Get the state i ---*/

      VelMagnitude2_i = 0.0; ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
        ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
        VelMagnitude2_i += Velocity_i[iDim]*Velocity_i[iDim];
      }
      Density_i = node[iPoint]->GetDensity();
      Energy_i = node[iPoint]->GetEnergy();

      /*--- Compute the boundary state b ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        Velocity_b[iDim] = Velocity_i[iDim] - ProjVelocity_i * UnitNormal[iDim]; //Force the velocity to be tangential to the surface.

      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
        for (iDim = 0; iDim < nDim; iDim++) Velocity_b[iDim] += GridVel[iDim] - ProjGridVel * UnitNormal[iDim];
      }

      VelMagnitude2_b = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        VelMagnitude2_b += Velocity_b[iDim] * Velocity_b[iDim];

      /*--- Compute the residual ---*/

      turb_ke = 0.0;
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);

      Density_b = Density_i;
      StaticEnergy_b = Energy_i - 0.5 * VelMagnitude2_i - turb_ke;
      Energy_b = StaticEnergy_b + 0.5 * VelMagnitude2_b + turb_ke;

      FluidModel->SetTDState_rhoe(Density_b, StaticEnergy_b);
      Kappa_b = FluidModel->GetdPde_rho() / Density_b;
      Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
      Pressure_b = FluidModel->GetPressure();
      Enthalpy_b = Energy_b + Pressure_b/Density_b;

      numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, NormalArea, Residual);

      /*--- Grid velocity correction to the energy term ---*/
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
        Residual[nVar-1] += Pressure_b*ProjGridVel*Area;
      }

      /*--- Add the Reynolds stress tensor contribution ---*/

      if (tkeNeeded) {
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[iDim+1] += (2.0/3.0)*Density_b*turb_ke*NormalArea[iDim];
      }
      
      /*--- Add value to the residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Form Jacobians for implicit computations ---*/
      
      if (implicit) {
        
        /*--- Initialize Jacobian ---*/
        
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
        
        /*--- Compute DubDu ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++)
            DubDu[iVar][jVar]= 0.0;
          DubDu[iVar][iVar]= 1.0;
        }

        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim<nDim; jDim++)
            DubDu[iDim+1][jDim+1] -= UnitNormal[iDim]*UnitNormal[jDim];
        DubDu[nVar-1][0] += 0.5*ProjVelocity_i*ProjVelocity_i;
        for (iDim = 0; iDim < nDim; iDim++) {
          DubDu[nVar-1][iDim+1] -= ProjVelocity_i*UnitNormal[iDim];
        }

        /*--- Compute flux Jacobian in state b ---*/

        numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, NormalArea, 1, Jacobian_b);

        // Check for grid movement, should be already considered since Jacobian b is computed from u_b
        // if (grid_movement) {
        // Jacobian_b[nVar-1][0] += 0.5*ProjGridVel*ProjGridVel;
        // for (iDim = 0; iDim < nDim; iDim++)
        // Jacobian_b[nVar-1][iDim+1] -= ProjGridVel * UnitNormal[iDim];
        // }

        /*--- Compute numerical flux Jacobian at node i ---*/

        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            for (kVar = 0; kVar < nVar; kVar++)
              Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];

        /*--- Add the Jacobian to the sparse matrix ---*/

        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      }
    }
  }
  
  delete [] Normal;
  delete [] NormalArea;
  delete [] Velocity_b;
  delete [] Velocity_i;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_b[iVar];
    delete [] DubDu[iVar];
  }
  delete [] Jacobian_b;
  delete [] DubDu;
  
}

void CEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  
  su2double *GridVel;
  su2double Area, UnitNormal[3] = {0.0,0.0,0.0};
  su2double Density, Pressure, Energy,  Velocity[3] = {0.0,0.0,0.0};
  su2double Density_Bound, Pressure_Bound, Vel_Bound[3] = {0.0,0.0,0.0};
  su2double Density_Infty, Pressure_Infty, Vel_Infty[3] = {0.0,0.0,0.0};
  su2double SoundSpeed, Entropy, Velocity2, Vn;
  su2double SoundSpeed_Bound, Entropy_Bound, Vel2_Bound, Vn_Bound;
  su2double SoundSpeed_Infty, Entropy_Infty, Vel2_Infty, Vn_Infty, Qn_Infty;
  su2double RiemannPlus, RiemannMinus;
  su2double *V_infty, *V_domain;
  
  su2double Gas_Constant     = config->GetGas_ConstantND();
  
  bool implicit       = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  bool grid_movement  = config->GetGrid_Movement();
  bool viscous        = config->GetViscous();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS ) ||
                     (config->GetKind_Solver() == DISC_ADJ_RANS))
                    && (config->GetKind_Turb_Model() == SST));
    
  su2double *Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Allocate the value at the infinity ---*/
    V_infty = GetCharacPrimVar(val_marker, iVertex);
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();

      /*--- Construct solution state at infinity for compressible flow by
         using Riemann invariants, and then impose a weak boundary condition
         by computing the flux using this new state for U. See CFD texts by
         Hirsch or Blazek for more detail. Adapted from an original
         implementation in the Stanford University multi-block (SUmb) solver
         in the routine bcFarfield.f90 written by Edwin van der Weide,
         last modified 06-12-2005. First, compute the unit normal at the
         boundary nodes. ---*/

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Store primitive variables (density, velocities, velocity squared,
         energy, pressure, and sound speed) at the boundary node, and set some
         other quantities for clarity. Project the current flow velocity vector
         at this boundary node into the local normal direction, i.e. compute
         v_bound.n.  ---*/

      Density_Bound = V_domain[nDim+2];
      Vel2_Bound = 0.0; Vn_Bound = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Bound[iDim] = V_domain[iDim+1];
        Vel2_Bound     += Vel_Bound[iDim]*Vel_Bound[iDim];
        Vn_Bound       += Vel_Bound[iDim]*UnitNormal[iDim];
      }
      Pressure_Bound   = node[iPoint]->GetPressure();
      SoundSpeed_Bound = sqrt(Gamma*Pressure_Bound/Density_Bound);
      Entropy_Bound    = pow(Density_Bound, Gamma)/Pressure_Bound;

      /*--- Store the primitive variable state for the freestream. Project
         the freestream velocity vector into the local normal direction,
         i.e. compute v_infty.n. ---*/

      Density_Infty = GetDensity_Inf();
      Vel2_Infty = 0.0; Vn_Infty = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Infty[iDim] = GetVelocity_Inf(iDim);
        Vel2_Infty     += Vel_Infty[iDim]*Vel_Infty[iDim];
        Vn_Infty       += Vel_Infty[iDim]*UnitNormal[iDim];
      }
      Pressure_Infty   = GetPressure_Inf();
      SoundSpeed_Infty = sqrt(Gamma*Pressure_Infty/Density_Infty);
      Entropy_Infty    = pow(Density_Infty, Gamma)/Pressure_Infty;

      /*--- Adjust the normal freestream velocity for grid movement ---*/

      Qn_Infty = Vn_Infty;
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          Qn_Infty -= GridVel[iDim]*UnitNormal[iDim];
      }

      /*--- Compute acoustic Riemann invariants: R = u.n +/- 2c/(gamma-1).
         These correspond with the eigenvalues (u+c) and (u-c), respectively,
         which represent the acoustic waves. Positive characteristics are
         incoming, and a physical boundary condition is imposed (freestream
         state). This occurs when either (u.n+c) > 0 or (u.n-c) > 0. Negative
         characteristics are leaving the domain, and numerical boundary
         conditions are required by extrapolating from the interior state
         using the Riemann invariants. This occurs when (u.n+c) < 0 or
         (u.n-c) < 0. Note that grid movement is taken into account when
         checking the sign of the eigenvalue. ---*/

      /*--- Check whether (u.n+c) is greater or less than zero ---*/

      if (Qn_Infty > -SoundSpeed_Infty) {
        /*--- Subsonic inflow or outflow ---*/
        RiemannPlus = Vn_Bound + 2.0*SoundSpeed_Bound/Gamma_Minus_One;
      } else {
        /*--- Supersonic inflow ---*/
        RiemannPlus = Vn_Infty + 2.0*SoundSpeed_Infty/Gamma_Minus_One;
      }

      /*--- Check whether (u.n-c) is greater or less than zero ---*/

      if (Qn_Infty > SoundSpeed_Infty) {
        /*--- Supersonic outflow ---*/
        RiemannMinus = Vn_Bound - 2.0*SoundSpeed_Bound/Gamma_Minus_One;
      } else {
        /*--- Subsonic outflow ---*/
        RiemannMinus = Vn_Infty - 2.0*SoundSpeed_Infty/Gamma_Minus_One;
      }

      /*--- Compute a new value for the local normal velocity and speed of
         sound from the Riemann invariants. ---*/

      Vn = 0.5 * (RiemannPlus + RiemannMinus);
      SoundSpeed = 0.25 * (RiemannPlus - RiemannMinus)*Gamma_Minus_One;

      /*--- Construct the primitive variable state at the boundary for
         computing the flux for the weak boundary condition. The values
         that we choose to construct the solution (boundary or freestream)
         depend on whether we are at an inflow or outflow. At an outflow, we
         choose boundary information (at most one characteristic is incoming),
         while at an inflow, we choose infinity values (at most one
         characteristic is outgoing). ---*/

      if (Qn_Infty > 0.0)   {
        /*--- Outflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Bound[iDim] + (Vn-Vn_Bound)*UnitNormal[iDim];
        Entropy = Entropy_Bound;
      } else  {
        /*--- Inflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Infty[iDim] + (Vn-Vn_Infty)*UnitNormal[iDim];
        Entropy = Entropy_Infty;
      }

      /*--- Recompute the primitive variables. ---*/

      Density = pow(Entropy*SoundSpeed*SoundSpeed/Gamma,1.0/Gamma_Minus_One);
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Pressure = Density*SoundSpeed*SoundSpeed/Gamma;
      Energy   = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;
      if (tkeNeeded) Energy += GetTke_Inf();

      /*--- Store new primitive state for computing the flux. ---*/

      V_infty[0] = Pressure/(Gas_Constant*Density);
      for (iDim = 0; iDim < nDim; iDim++)
        V_infty[iDim+1] = Velocity[iDim];
      V_infty[nDim+1] = Pressure;
      V_infty[nDim+2] = Density;
      V_infty[nDim+3] = Energy + Pressure/Density;


      
      /*--- Set various quantities in the numerics class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      if (grid_movement) {
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      }
      
      /*--- Compute the convective residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Convective Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous residual contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_infty[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_infty[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_infty);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),
                                          node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                              solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update viscous residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Viscous Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
}

void CEulerSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container,
                              CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iVar, jVar, kVar;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double P_Total, T_Total, P_static, T_static, Rho_static, *Mach, *Flow_Dir, Area, UnitNormal[3];
  su2double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, StaticEnergy_b, Density_b, Kappa_b, Chi_b, Pressure_b, Temperature_b;
  su2double *Velocity_e, Velocity2_e, VelMag_e, Enthalpy_e, Entropy_e, Energy_e = 0.0, StaticEnthalpy_e, StaticEnergy_e, Density_e = 0.0, Pressure_e;
  su2double *Velocity_i, Velocity2_i, Enthalpy_i, Energy_i, StaticEnergy_i, Density_i, Kappa_i, Chi_i, Pressure_i, SoundSpeed_i;
  su2double ProjVelocity_i;
  su2double **P_Tensor, **invP_Tensor, *Lambda_i, **Jacobian_b, **DubDu, *dw, *u_e, *u_i, *u_b;
  su2double *gridVel;
  su2double *V_boundary, *V_domain, *S_boundary, *S_domain;
  
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  
  su2double *Normal, *FlowDirMix, TangVelocity, NormalVelocity;
  Normal = new su2double[nDim];

  Velocity_i = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  Velocity_e = new su2double[nDim];
  FlowDirMix = new su2double[nDim];
  Lambda_i = new su2double[nVar];
  u_i = new su2double[nVar];
  u_e = new su2double[nVar];
  u_b = new su2double[nVar];
  dw = new su2double[nVar];
  
  S_boundary = new su2double[8];
  
  P_Tensor = new su2double*[nVar];
  invP_Tensor = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
  {
    P_Tensor[iVar] = new su2double[nVar];
    invP_Tensor[iVar] = new su2double[nVar];
  }
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    V_boundary= GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Compute the internal state u_i ---*/
      Velocity2_i = 0;
      for (iDim=0; iDim < nDim; iDim++)
      {
        Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
        Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
      }
      
      
      Density_i = node[iPoint]->GetDensity();
      
      Energy_i = node[iPoint]->GetEnergy();
      StaticEnergy_i = Energy_i - 0.5*Velocity2_i;
      
      FluidModel->SetTDState_rhoe(Density_i, StaticEnergy_i);
      
      Pressure_i = FluidModel->GetPressure();
      Enthalpy_i = Energy_i + Pressure_i/Density_i;
      
      SoundSpeed_i = FluidModel->GetSoundSpeed();
      
      Kappa_i = FluidModel->GetdPde_rho() / Density_i;
      Chi_i = FluidModel->GetdPdrho_e() - Kappa_i * StaticEnergy_i;
      
      ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      
      /*--- Build the external state u_e from boundary data and internal node ---*/
      
      switch(config->GetKind_Data_Riemann(Marker_Tag))
      {
      case TOTAL_CONDITIONS_PT:

        /*--- Retrieve the specified total conditions for this boundary. ---*/
        if (gravity) P_Total = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;/// check in which case is true (only freesurface?)
        else P_Total  = config->GetRiemann_Var1(Marker_Tag);
        T_Total  = config->GetRiemann_Var2(Marker_Tag);
        Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        P_Total /= config->GetPressure_Ref();
        T_Total /= config->GetTemperature_Ref();

        /* --- Computes the total state --- */
        FluidModel->SetTDState_PT(P_Total, T_Total);
        Enthalpy_e = FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
        Entropy_e = FluidModel->GetEntropy();

        /* --- Compute the boundary state u_e --- */
        Velocity2_e = Velocity2_i;
        if (nDim == 2){
          NormalVelocity= -sqrt(Velocity2_e)*Flow_Dir[0];
          TangVelocity= -sqrt(Velocity2_e)*Flow_Dir[1];
          Velocity_e[0]= UnitNormal[0]*NormalVelocity - UnitNormal[1]*TangVelocity;
          Velocity_e[1]= UnitNormal[1]*NormalVelocity + UnitNormal[0]*TangVelocity;
        }else{
          for (iDim = 0; iDim < nDim; iDim++)
            Velocity_e[iDim] = sqrt(Velocity2_e)*Flow_Dir[iDim];
        }
        StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
        FluidModel->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
        Density_e = FluidModel->GetDensity();
        StaticEnergy_e = FluidModel->GetStaticEnergy();
        Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
        if (tkeNeeded) Energy_e += GetTke_Inf();
        break;

      case STATIC_SUPERSONIC_INFLOW_PT:

        /*--- Retrieve the specified total conditions for this boundary. ---*/
        if (gravity) P_static = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;/// check in which case is true (only freesurface?)
        else P_static  = config->GetRiemann_Var1(Marker_Tag);
        T_static  = config->GetRiemann_Var2(Marker_Tag);
        Mach = config->GetRiemann_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        P_static /= config->GetPressure_Ref();
        T_static /= config->GetTemperature_Ref();

        /* --- Computes the total state --- */
        FluidModel->SetTDState_PT(P_static, T_static);

        /* --- Compute the boundary state u_e --- */
        Velocity2_e = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity_e[iDim] = Mach[iDim]*FluidModel->GetSoundSpeed();
          Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
        }
        Density_e = FluidModel->GetDensity();
        StaticEnergy_e = FluidModel->GetStaticEnergy();
        Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
        if (tkeNeeded) Energy_e += GetTke_Inf();
        break;

      case STATIC_SUPERSONIC_INFLOW_PD:

        /*--- Retrieve the specified total conditions for this boundary. ---*/

        if (gravity) P_static = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;/// check in which case is true (only freesurface?)
        else P_static  = config->GetRiemann_Var1(Marker_Tag);
        Rho_static  = config->GetRiemann_Var2(Marker_Tag);
        Mach = config->GetRiemann_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        P_static /= config->GetPressure_Ref();
        Rho_static /= config->GetDensity_Ref();

        /* --- Computes the total state --- */
        FluidModel->SetTDState_Prho(P_static, Rho_static);

        /* --- Compute the boundary state u_e --- */
        Velocity2_e = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity_e[iDim] = Mach[iDim]*FluidModel->GetSoundSpeed();
          Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
        }
        Density_e = FluidModel->GetDensity();
        StaticEnergy_e = FluidModel->GetStaticEnergy();
        Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
        if (tkeNeeded) Energy_e += GetTke_Inf();
        break;

      case DENSITY_VELOCITY:

        /*--- Retrieve the specified density and velocity magnitude ---*/
        Density_e  = config->GetRiemann_Var1(Marker_Tag);
        VelMag_e   = config->GetRiemann_Var2(Marker_Tag);
        Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        Density_e /= config->GetDensity_Ref();
        VelMag_e /= config->GetVelocity_Ref();

        for (iDim = 0; iDim < nDim; iDim++)
          Velocity_e[iDim] = VelMag_e*Flow_Dir[iDim];
        Energy_e = Energy_i;
        break;

      case STATIC_PRESSURE:

        /*--- Retrieve the staic pressure for this boundary. ---*/
        Pressure_e = config->GetRiemann_Var1(Marker_Tag);
        Pressure_e /= config->GetPressure_Ref();
        Density_e = Density_i;

        /* --- Compute the boundary state u_e --- */
        FluidModel->SetTDState_Prho(Pressure_e, Density_e);
        Velocity2_e = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity_e[iDim] = Velocity_i[iDim];
          Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
        }
        Energy_e = FluidModel->GetStaticEnergy() + 0.5*Velocity2_e;
        break;

      default:
        SU2_MPI::Error("Invalid Riemann input!", CURRENT_FUNCTION);
        break;
      }
      
      /*--- Compute P (matrix of right eigenvectors) ---*/
      conv_numerics->GetPMatrix(&Density_i, Velocity_i, &SoundSpeed_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, P_Tensor);
      
      /*--- Compute inverse P (matrix of left eigenvectors)---*/
      conv_numerics->GetPMatrix_inv(invP_Tensor, &Density_i, Velocity_i, &SoundSpeed_i, &Chi_i, &Kappa_i, UnitNormal);
      
      /*--- eigenvalues contribution due to grid motion ---*/
      if (grid_movement) {
        gridVel = geometry->node[iPoint]->GetGridVel();
        
        su2double ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel   += gridVel[iDim]*UnitNormal[iDim];
        ProjVelocity_i -= ProjGridVel;
      }
      
      /*--- Flow eigenvalues ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Lambda_i[iDim] = ProjVelocity_i;
      Lambda_i[nVar-2] = ProjVelocity_i + SoundSpeed_i;
      Lambda_i[nVar-1] = ProjVelocity_i - SoundSpeed_i;
      
      /*--- Compute the boundary state u_e ---*/
      u_e[0] = Density_e;
      for (iDim = 0; iDim < nDim; iDim++)
        u_e[iDim+1] = Velocity_e[iDim]*Density_e;
      u_e[nVar-1] = Energy_e*Density_e;
      
      /*--- Compute the boundary state u_i ---*/
      u_i[0] = Density_i;
      for (iDim = 0; iDim < nDim; iDim++)
        u_i[iDim+1] = Velocity_i[iDim]*Density_i;
      u_i[nVar-1] = Energy_i*Density_i;
      
      /*--- Compute the characteristic jumps ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        dw[iVar] = 0;
        for (jVar = 0; jVar < nVar; jVar++)
          dw[iVar] += invP_Tensor[iVar][jVar] * (u_e[jVar] - u_i[jVar]);
        
      }
      
      /*--- Compute the boundary state u_b using characteristics ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        u_b[iVar] = u_i[iVar];
        
        for (jVar = 0; jVar < nVar; jVar++)
        {
          if (Lambda_i[jVar] < 0)
          {
            u_b[iVar] += P_Tensor[iVar][jVar]*dw[jVar];
            
          }
        }
      }
      
      
      /*--- Compute the thermodynamic state in u_b ---*/
      Density_b = u_b[0];
      Velocity2_b = 0;
      for (iDim = 0; iDim < nDim; iDim++)
      {
        Velocity_b[iDim] = u_b[iDim+1]/Density_b;
        Velocity2_b += Velocity_b[iDim]*Velocity_b[iDim];
      }
      Energy_b = u_b[nVar-1]/Density_b;
      StaticEnergy_b = Energy_b - 0.5*Velocity2_b;
      FluidModel->SetTDState_rhoe(Density_b, StaticEnergy_b);
      Pressure_b = FluidModel->GetPressure();
      Temperature_b = FluidModel->GetTemperature();
      Enthalpy_b = Energy_b + Pressure_b/Density_b;
      Kappa_b = FluidModel->GetdPde_rho() / Density_b;
      Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
      
      /*--- Compute the residuals ---*/
      conv_numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, Normal, Residual);
      
      /*--- Residual contribution due to grid motion ---*/
      if (grid_movement) {
        gridVel = geometry->node[iPoint]->GetGridVel();
        su2double projVelocity = 0.0;
        
        for (iDim = 0; iDim < nDim; iDim++)
          projVelocity +=  gridVel[iDim]*Normal[iDim];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] -= projVelocity *(u_b[iVar]);
      }
      
      if (implicit) {
        
        Jacobian_b = new su2double*[nVar];
        DubDu = new su2double*[nVar];
        for (iVar = 0; iVar < nVar; iVar++)
        {
          Jacobian_b[iVar] = new su2double[nVar];
          DubDu[iVar] = new su2double[nVar];
        }
        
        /*--- Initialize DubDu to unit matrix---*/
        
        for (iVar = 0; iVar < nVar; iVar++)
        {
          for (jVar = 0; jVar < nVar; jVar++)
            DubDu[iVar][jVar]= 0;
          
          DubDu[iVar][iVar]= 1;
        }
        
        /*--- Compute DubDu -= RNL---*/
        for (iVar=0; iVar<nVar; iVar++)
        {
          for (jVar=0; jVar<nVar; jVar++)
          {
            for (kVar=0; kVar<nVar; kVar++)
            {
              if (Lambda_i[kVar]<0)
                DubDu[iVar][jVar] -= P_Tensor[iVar][kVar] * invP_Tensor[kVar][jVar];
            }
          }
        }
        
        /*--- Compute flux Jacobian in state b ---*/
        conv_numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, Normal, 1.0, Jacobian_b);
        
        /*--- Jacobian contribution due to grid motion ---*/
        if (grid_movement)
        {
          gridVel = geometry->node[iPoint]->GetGridVel();
          su2double projVelocity = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            projVelocity +=  gridVel[iDim]*Normal[iDim];
          for (iVar = 0; iVar < nVar; iVar++) {
            Residual[iVar] -= projVelocity *(u_b[iVar]);
            Jacobian_b[iVar][iVar] -= projVelocity;
          }
          
        }
        
        /*--- initiate Jacobian_i to zero matrix ---*/
        for (iVar=0; iVar<nVar; iVar++)
          for (jVar=0; jVar<nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        
        /*--- Compute numerical flux Jacobian at node i ---*/
        for (iVar=0; iVar<nVar; iVar++) {
          for (jVar=0; jVar<nVar; jVar++) {
            for (kVar=0; kVar<nVar; kVar++) {
              Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];
            }
          }
        }
        
        for (iVar = 0; iVar < nVar; iVar++) {
          delete [] Jacobian_b[iVar];
          delete [] DubDu[iVar];
        }
        delete [] Jacobian_b;
        delete [] DubDu;
      }
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      if (viscous) {
        
        /*--- Primitive variables, using the derived quantities ---*/
        V_boundary[0] = Temperature_b;
        for (iDim = 0; iDim < nDim; iDim++)
          V_boundary[iDim+1] = Velocity_b[iDim];
        V_boundary[nDim+1] = Pressure_b;
        V_boundary[nDim+2] = Density_b;
        V_boundary[nDim+3] = Enthalpy_b;
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_boundary[nDim+5] = FluidModel->GetLaminarViscosity();
        V_boundary[nDim+6] = node[iPoint]->GetEddyViscosity();
        V_boundary[nDim+7] = FluidModel->GetThermalConductivity();
        V_boundary[nDim+8] = FluidModel->GetCp();
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_boundary);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Secondary variables ---*/
        S_domain = node[iPoint]->GetSecondary();
        
        /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/
        
        S_boundary[0]= FluidModel->GetdPdrho_e();
        S_boundary[1]= FluidModel->GetdPde_rho();
        
        S_boundary[2]= FluidModel->GetdTdrho_e();
        S_boundary[3]= FluidModel->GetdTde_rho();
        
        /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/
        
        S_boundary[4]= FluidModel->Getdmudrho_T();
        S_boundary[5]= FluidModel->GetdmudT_rho();
        
        S_boundary[6]= FluidModel->Getdktdrho_T();
        S_boundary[7]= FluidModel->GetdktdT_rho();
        
        visc_numerics->SetSecondary(S_domain, S_boundary);
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  delete [] Velocity_e;
  delete [] Velocity_b;
  delete [] Velocity_i;
  delete [] FlowDirMix;
  
  delete [] S_boundary;
  delete [] Lambda_i;
  delete [] u_i;
  delete [] u_e;
  delete [] u_b;
  delete [] dw;
  
  
  for (iVar = 0; iVar < nVar; iVar++)
  {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}


void CEulerSolver::BC_TurboRiemann(CGeometry *geometry, CSolver **solver_container,
    CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iVar, jVar, kVar, iSpan;
  unsigned long iPoint, Point_Normal, oldVertex;
  long iVertex;
  su2double P_Total, T_Total, *Flow_Dir;
  su2double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, StaticEnergy_b, Density_b, Kappa_b, Chi_b, Pressure_b, Temperature_b;
  su2double *Velocity_e, Velocity2_e, Enthalpy_e, Entropy_e, Energy_e = 0.0, StaticEnthalpy_e, StaticEnergy_e, Density_e = 0.0, Pressure_e;
  su2double *Velocity_i, Velocity2_i, Enthalpy_i, Energy_i, StaticEnergy_i, Density_i, Kappa_i, Chi_i, Pressure_i, SoundSpeed_i;
  su2double ProjVelocity_i;
  su2double **P_Tensor, **invP_Tensor, *Lambda_i, **Jacobian_b, **DubDu, *dw, *u_e, *u_i, *u_b;
  su2double *gridVel;
  su2double *V_boundary, *V_domain, *S_boundary, *S_domain;
  su2double AverageEnthalpy, AverageEntropy;
  unsigned short  iZone     = config->GetiZone();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  unsigned short nSpanWiseSections = geometry->GetnSpanWiseSections(config->GetMarker_All_TurbomachineryFlag(val_marker));
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
      (config->GetKind_Turb_Model() == SST));

  su2double *Normal, *turboNormal, *UnitNormal, *FlowDirMix, FlowDirMixMag, *turboVelocity;
  Normal = new su2double[nDim];
  turboNormal 	= new su2double[nDim];
  UnitNormal 	= new su2double[nDim];

  Velocity_i = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  Velocity_e = new su2double[nDim];
  turboVelocity = new su2double[nDim];
  FlowDirMix = new su2double[nDim];
  Lambda_i = new su2double[nVar];
  u_i = new su2double[nVar];
  u_e = new su2double[nVar];
  u_b = new su2double[nVar];
  dw = new su2double[nVar];

  S_boundary = new su2double[8];

  P_Tensor = new su2double*[nVar];
  invP_Tensor = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
  {
    P_Tensor[iVar] = new su2double[nVar];
    invP_Tensor[iVar] = new su2double[nVar];
  }

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iSpan= 0; iSpan < nSpanWiseSections; iSpan++){
    for (iVertex = 0; iVertex < geometry->nVertexSpan[val_marker][iSpan]; iVertex++) {

      /*--- using the other vertex information for retrieving some information ---*/
      oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();
      V_boundary= GetCharacPrimVar(val_marker, oldVertex);

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention),
       * 		this normal is scaled with the area of the face of the element  ---*/
      geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- find the node related to the vertex ---*/
      iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();
      if (geometry->node[iPoint]->GetDomain()){

        /*--- Normalize Normal vector for this vertex (already for outward convention) ---*/
        geometry->turbovertex[val_marker][iSpan][iVertex]->GetNormal(UnitNormal);
        geometry->turbovertex[val_marker][iSpan][iVertex]->GetTurboNormal(turboNormal);

        /*--- Retrieve solution at this boundary node ---*/
        V_domain = node[iPoint]->GetPrimitive();

        /* --- Compute the internal state u_i --- */
        Velocity2_i = 0;
        for (iDim=0; iDim < nDim; iDim++)
        {
          Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
          Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
        }

        Density_i = node[iPoint]->GetDensity();

        Energy_i = node[iPoint]->GetEnergy();
        StaticEnergy_i = Energy_i - 0.5*Velocity2_i;

        FluidModel->SetTDState_rhoe(Density_i, StaticEnergy_i);

        Pressure_i = FluidModel->GetPressure();
        Enthalpy_i = Energy_i + Pressure_i/Density_i;

        SoundSpeed_i = FluidModel->GetSoundSpeed();

        Kappa_i = FluidModel->GetdPde_rho() / Density_i;
        Chi_i = FluidModel->GetdPdrho_e() - Kappa_i * StaticEnergy_i;

        ProjVelocity_i = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];

        /*--- Build the external state u_e from boundary data and internal node ---*/

        switch(config->GetKind_Data_Riemann(Marker_Tag))
        {
          //TODO(turbo), generilize for 3D case
          //TODO(turbo), generilize for Inlet and Outlet in for backflow treatment
          //TODO(turbo), implement not uniform inlet and radial equilibrium for the outlet
          case TOTAL_CONDITIONS_PT:

            /*--- Retrieve the specified total conditions for this boundary. ---*/
            if (gravity) P_Total = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;/// check in which case is true (only freesurface?)
            else P_Total  = config->GetRiemann_Var1(Marker_Tag);
            T_Total  = config->GetRiemann_Var2(Marker_Tag);
            Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);

            /*--- Non-dim. the inputs if necessary. ---*/
            P_Total /= config->GetPressure_Ref();
            T_Total /= config->GetTemperature_Ref();

            /* --- Computes the total state --- */
            FluidModel->SetTDState_PT(P_Total, T_Total);
            Enthalpy_e = FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
            Entropy_e = FluidModel->GetEntropy();

            /* --- Compute the boundary state u_e --- */
            Velocity2_e = Velocity2_i;
            for (iDim = 0; iDim < nDim; iDim++)
              turboVelocity[iDim] = sqrt(Velocity2_e)*Flow_Dir[iDim];
            ComputeBackVelocity(turboVelocity,turboNormal, Velocity_e, config->GetMarker_All_TurbomachineryFlag(val_marker),config->GetKind_TurboMachinery(iZone));
            StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
            FluidModel->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
            Density_e = FluidModel->GetDensity();
            StaticEnergy_e = FluidModel->GetStaticEnergy();
            Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
            if (tkeNeeded) Energy_e += GetTke_Inf();
            break;

          case MIXING_IN:

            /* --- compute total averaged quantities ---*/
            FluidModel->SetTDState_Prho(ExtAveragePressure[val_marker][iSpan], ExtAverageDensity[val_marker][iSpan]);
            AverageEnthalpy = FluidModel->GetStaticEnergy() + ExtAveragePressure[val_marker][iSpan]/ExtAverageDensity[val_marker][iSpan];
            AverageEntropy  = FluidModel->GetEntropy();

            FlowDirMixMag = 0;
            for (iDim = 0; iDim < nDim; iDim++)
              FlowDirMixMag += ExtAverageTurboVelocity[val_marker][iSpan][iDim]*ExtAverageTurboVelocity[val_marker][iSpan][iDim];
            for (iDim = 0; iDim < nDim; iDim++){
              FlowDirMix[iDim] = ExtAverageTurboVelocity[val_marker][iSpan][iDim]/sqrt(FlowDirMixMag);
            }


            /* --- Computes the total state --- */
            Enthalpy_e = AverageEnthalpy;
            Entropy_e = AverageEntropy;

            /* --- Compute the boundary state u_e --- */
            Velocity2_e = Velocity2_i;
            for (iDim = 0; iDim < nDim; iDim++){
              turboVelocity[iDim] = sqrt(Velocity2_e)*FlowDirMix[iDim];

            }
            ComputeBackVelocity(turboVelocity,turboNormal, Velocity_e, config->GetMarker_All_TurbomachineryFlag(val_marker),config->GetKind_TurboMachinery(iZone));

            StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
            FluidModel->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
            Density_e = FluidModel->GetDensity();
            StaticEnergy_e = FluidModel->GetStaticEnergy();
            Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
            //				if (tkeNeeded) Energy_e += GetTke_Inf();
            break;


          case MIXING_OUT:

            /*--- Retrieve the staic pressure for this boundary. ---*/
            Pressure_e = ExtAveragePressure[val_marker][iSpan];
            Density_e = Density_i;

            /* --- Compute the boundary state u_e --- */
            FluidModel->SetTDState_Prho(Pressure_e, Density_e);
            Velocity2_e = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Velocity_e[iDim] = Velocity_i[iDim];
              Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
            }
            Energy_e = FluidModel->GetStaticEnergy() + 0.5*Velocity2_e;
            break;

          case STATIC_PRESSURE:

            /*--- Retrieve the staic pressure for this boundary. ---*/
            Pressure_e = config->GetRiemann_Var1(Marker_Tag);
            Pressure_e /= config->GetPressure_Ref();
            Density_e = Density_i;

            /* --- Compute the boundary state u_e --- */
            FluidModel->SetTDState_Prho(Pressure_e, Density_e);
            Velocity2_e = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Velocity_e[iDim] = Velocity_i[iDim];
              Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
            }
            Energy_e = FluidModel->GetStaticEnergy() + 0.5*Velocity2_e;
            break;


          case RADIAL_EQUILIBRIUM:

            /*--- Retrieve the staic pressure for this boundary. ---*/
            Pressure_e = RadialEquilibriumPressure[val_marker][iSpan];
            Density_e = Density_i;

            /* --- Compute the boundary state u_e --- */
            FluidModel->SetTDState_Prho(Pressure_e, Density_e);
            Velocity2_e = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Velocity_e[iDim] = Velocity_i[iDim];
              Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
            }
            Energy_e = FluidModel->GetStaticEnergy() + 0.5*Velocity2_e;
            break;

          default:
            SU2_MPI::Error("Invalid Riemann input!", CURRENT_FUNCTION);
            break;
        }

        /*--- Compute P (matrix of right eigenvectors) ---*/
        conv_numerics->GetPMatrix(&Density_i, Velocity_i, &SoundSpeed_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, P_Tensor);

        /*--- Compute inverse P (matrix of left eigenvectors)---*/
        conv_numerics->GetPMatrix_inv(invP_Tensor, &Density_i, Velocity_i, &SoundSpeed_i, &Chi_i, &Kappa_i, UnitNormal);

        /*--- eigenvalues contribution due to grid motion ---*/
        if (grid_movement){
          gridVel = geometry->node[iPoint]->GetGridVel();

          su2double ProjGridVel = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel   += gridVel[iDim]*UnitNormal[iDim];
          ProjVelocity_i -= ProjGridVel;
        }

        /*--- Flow eigenvalues ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Lambda_i[iDim] = ProjVelocity_i;
        Lambda_i[nVar-2] = ProjVelocity_i + SoundSpeed_i;
        Lambda_i[nVar-1] = ProjVelocity_i - SoundSpeed_i;

        /* --- Compute the boundary state u_e --- */
        u_e[0] = Density_e;
        for (iDim = 0; iDim < nDim; iDim++)
          u_e[iDim+1] = Velocity_e[iDim]*Density_e;
        u_e[nVar-1] = Energy_e*Density_e;

        /* --- Compute the boundary state u_i --- */
        u_i[0] = Density_i;
        for (iDim = 0; iDim < nDim; iDim++)
          u_i[iDim+1] = Velocity_i[iDim]*Density_i;
        u_i[nVar-1] = Energy_i*Density_i;

        /*--- Compute the characteristic jumps ---*/
        for (iVar = 0; iVar < nVar; iVar++)
        {
          dw[iVar] = 0;
          for (jVar = 0; jVar < nVar; jVar++)
            dw[iVar] += invP_Tensor[iVar][jVar] * (u_e[jVar] - u_i[jVar]);

        }

        /*--- Compute the boundary state u_b using characteristics ---*/
        for (iVar = 0; iVar < nVar; iVar++)
        {
          u_b[iVar] = u_i[iVar];

          for (jVar = 0; jVar < nVar; jVar++)
          {
            if (Lambda_i[jVar] < 0)
            {
              u_b[iVar] += P_Tensor[iVar][jVar]*dw[jVar];

            }
          }
        }


        /*--- Compute the thermodynamic state in u_b ---*/
        Density_b = u_b[0];
        Velocity2_b = 0;
        for (iDim = 0; iDim < nDim; iDim++)
        {
          Velocity_b[iDim] = u_b[iDim+1]/Density_b;
          Velocity2_b += Velocity_b[iDim]*Velocity_b[iDim];
        }
        Energy_b = u_b[nVar-1]/Density_b;
        StaticEnergy_b = Energy_b - 0.5*Velocity2_b;
        FluidModel->SetTDState_rhoe(Density_b, StaticEnergy_b);
        Pressure_b = FluidModel->GetPressure();
        Temperature_b = FluidModel->GetTemperature();
        Enthalpy_b = Energy_b + Pressure_b/Density_b;
        Kappa_b = FluidModel->GetdPde_rho() / Density_b;
        Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;

        /*--- Compute the residuals ---*/
        conv_numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, Normal, Residual);

        /*--- Residual contribution due to grid motion ---*/
        if (grid_movement) {
          gridVel = geometry->node[iPoint]->GetGridVel();
          su2double projVelocity = 0.0;

          for (iDim = 0; iDim < nDim; iDim++)
            projVelocity +=  gridVel[iDim]*Normal[iDim];
          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] -= projVelocity *(u_b[iVar]);
        }

        if (implicit) {

          Jacobian_b = new su2double*[nVar];
          DubDu = new su2double*[nVar];
          for (iVar = 0; iVar < nVar; iVar++)
          {
            Jacobian_b[iVar] = new su2double[nVar];
            DubDu[iVar] = new su2double[nVar];
          }

          /*--- Initialize DubDu to unit matrix---*/

          for (iVar = 0; iVar < nVar; iVar++)
          {
            for (jVar = 0; jVar < nVar; jVar++)
              DubDu[iVar][jVar]= 0;

            DubDu[iVar][iVar]= 1;
          }

          /*--- Compute DubDu -= RNL---*/
          for (iVar=0; iVar<nVar; iVar++)
          {
            for (jVar=0; jVar<nVar; jVar++)
            {
              for (kVar=0; kVar<nVar; kVar++)
              {
                if (Lambda_i[kVar]<0)
                  DubDu[iVar][jVar] -= P_Tensor[iVar][kVar] * invP_Tensor[kVar][jVar];
              }
            }
          }

          /*--- Compute flux Jacobian in state b ---*/
          conv_numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, Normal, 1.0, Jacobian_b);

          /*--- Jacobian contribution due to grid motion ---*/
          if (grid_movement)
          {
            gridVel = geometry->node[iPoint]->GetGridVel();
            su2double projVelocity = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              projVelocity +=  gridVel[iDim]*Normal[iDim];
            for (iVar = 0; iVar < nVar; iVar++){
              Residual[iVar] -= projVelocity *(u_b[iVar]);
              Jacobian_b[iVar][iVar] -= projVelocity;
            }

          }

          /*--- initiate Jacobian_i to zero matrix ---*/
          for (iVar=0; iVar<nVar; iVar++)
            for (jVar=0; jVar<nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;

          /*--- Compute numerical flux Jacobian at node i ---*/
          for (iVar=0; iVar<nVar; iVar++) {
            for (jVar=0; jVar<nVar; jVar++) {
              for (kVar=0; kVar<nVar; kVar++) {
                Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];
              }
            }
          }

          for (iVar = 0; iVar < nVar; iVar++) {
            delete [] Jacobian_b[iVar];
            delete [] DubDu[iVar];
          }
          delete [] Jacobian_b;
          delete [] DubDu;
        }

        /*--- Update residual value ---*/
        LinSysRes.AddBlock(iPoint, Residual);

        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

        /*--- Roe Turkel preconditioning, set the value of beta ---*/
        if (config->GetKind_Upwind() == TURKEL)
          node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());

        /*--- Viscous contribution ---*/
        if (viscous) {

          /*--- Primitive variables, using the derived quantities ---*/
          V_boundary[0] = Temperature_b;
          for (iDim = 0; iDim < nDim; iDim++)
            V_boundary[iDim+1] = Velocity_b[iDim];
          V_boundary[nDim+1] = Pressure_b;
          V_boundary[nDim+2] = Density_b;
          V_boundary[nDim+3] = Enthalpy_b;

          /*--- Set laminar and eddy viscosity at the infinity ---*/
          V_boundary[nDim+5] = FluidModel->GetLaminarViscosity();
          V_boundary[nDim+6] = node[iPoint]->GetEddyViscosity();
          V_boundary[nDim+7] = FluidModel->GetThermalConductivity();
          V_boundary[nDim+8] = FluidModel->GetCp();

          /*--- Set the normal vector and the coordinates ---*/
          visc_numerics->SetNormal(Normal);
          visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

          /*--- Primitive variables, and gradient ---*/
          visc_numerics->SetPrimitive(V_domain, V_boundary);
          visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

          /*--- Secondary variables ---*/
          S_domain = node[iPoint]->GetSecondary();

          /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/

          S_boundary[0]= FluidModel->GetdPdrho_e();
          S_boundary[1]= FluidModel->GetdPde_rho();

          S_boundary[2]= FluidModel->GetdTdrho_e();
          S_boundary[3]= FluidModel->GetdTde_rho();

          /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/

          S_boundary[4]= FluidModel->Getdmudrho_T();
          S_boundary[5]= FluidModel->GetdmudT_rho();

          S_boundary[6]= FluidModel->Getdktdrho_T();
          S_boundary[7]= FluidModel->GetdktdT_rho();

          visc_numerics->SetSecondary(S_domain, S_boundary);

          /*--- Turbulent kinetic energy ---*/
          if (config->GetKind_Turb_Model() == SST)
            visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));

          /*--- Compute and update residual ---*/
          visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
          LinSysRes.SubtractBlock(iPoint, Residual);

          /*--- Jacobian contribution for implicit integration ---*/
          if (implicit)
            Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

        }
      }
    }
}

  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  delete [] UnitNormal;
  delete [] turboNormal;
  delete [] turboVelocity;
  delete [] Velocity_e;
  delete [] Velocity_b;
  delete [] Velocity_i;
  delete [] FlowDirMix;

  delete [] S_boundary;
  delete [] Lambda_i;
  delete [] u_i;
  delete [] u_e;
  delete [] u_b;
  delete [] dw;


  for (iVar = 0; iVar < nVar; iVar++)
  {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;

}

void CEulerSolver::PreprocessBC_Giles(CGeometry *geometry, CConfig *config, CNumerics *conv_numerics, unsigned short marker_flag) {
  /* Implementation of Fuorier Transformations for non-regfelcting BC will come soon */
  su2double cj_inf,cj_out1, cj_out2, Density_i, Pressure_i, *turboNormal, *turboVelocity, *Velocity_i, AverageSoundSpeed;
  su2double *deltaprim, *cj, TwoPiThetaFreq_Pitch, pitch, theta, deltaTheta;
  unsigned short iMarker, iSpan, iMarkerTP, iDim;
  unsigned long  iPoint, kend_max, k;
  long iVertex, freq;
  unsigned short  iZone     = config->GetiZone();
  unsigned short nSpanWiseSections = geometry->GetnSpanWiseSections(marker_flag);
  turboNormal 	= new su2double[nDim];
  turboVelocity = new su2double[nDim];
  Velocity_i 		= new su2double[nDim];
  deltaprim     = new su2double[nVar];
  cj				    = new su2double[nVar];
  complex<su2double> I, cktemp_inf,cktemp_out1, cktemp_out2, expArg;
  I = complex<su2double>(0.0,1.0);

#ifdef HAVE_MPI
  su2double MyIm_inf, MyRe_inf, Im_inf, Re_inf, MyIm_out1, MyRe_out1, Im_out1, Re_out1, MyIm_out2, MyRe_out2, Im_out2, Re_out2;
#endif

  kend_max = geometry->GetnFreqSpanMax(marker_flag);
  for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){
    for(k=0; k < 2*kend_max+1; k++){
      freq = k - kend_max;
      cktemp_inf = complex<su2double>(0.0,0.0);
      cktemp_out1 = complex<su2double>(0.0,0.0);
      cktemp_out2 = complex<su2double>(0.0,0.0);
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
          if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
            if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
              for (iVertex = 0; iVertex < geometry->nVertexSpan[iMarker][iSpan]; iVertex++) {

                /*--- find the node related to the vertex ---*/
                iPoint = geometry->turbovertex[iMarker][iSpan][iVertex]->GetNode();

                geometry->turbovertex[iMarker][iSpan][iVertex]->GetTurboNormal(turboNormal);
                /*--- Compute the internal state _i ---*/

                Pressure_i = node[iPoint]->GetPressure();
                Density_i = node[iPoint]->GetDensity();
                for (iDim = 0; iDim < nDim; iDim++)
                {
                  Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
                }
                ComputeTurboVelocity(Velocity_i, turboNormal, turboVelocity, marker_flag, config->GetKind_TurboMachinery(iZone));

                if(nDim ==2){
                  deltaprim[0] = Density_i - AverageDensity[iMarker][iSpan];
                  deltaprim[1] = turboVelocity[0] - AverageTurboVelocity[iMarker][iSpan][0];
                  deltaprim[2] = turboVelocity[1] - AverageTurboVelocity[iMarker][iSpan][1];
                  deltaprim[3] = Pressure_i - AveragePressure[iMarker][iSpan];
                }
                else{
                  //Here 3d
                  deltaprim[0] = Density_i - AverageDensity[iMarker][iSpan];
                  deltaprim[1] = turboVelocity[0] - AverageTurboVelocity[iMarker][iSpan][0];
                  deltaprim[2] = turboVelocity[1] - AverageTurboVelocity[iMarker][iSpan][1];
                  deltaprim[3] = turboVelocity[2] - AverageTurboVelocity[iMarker][iSpan][2]; //New char
                  deltaprim[4] = Pressure_i - AveragePressure[iMarker][iSpan];
                }

                FluidModel->SetTDState_Prho(AveragePressure[iMarker][iSpan], AverageDensity[iMarker][iSpan]);
                AverageSoundSpeed = FluidModel->GetSoundSpeed();
                conv_numerics->GetCharJump(AverageSoundSpeed, AverageDensity[iMarker][iSpan], deltaprim, cj);

                /*-----this is only valid 2D ----*/
                if(nDim ==2){
                  cj_out1 = cj[1];
                  cj_out2 = cj[2];
                  cj_inf  = cj[3];
                }
                else{
                  //Here 3D
                  cj_out1 = cj[1];
                  cj_out2 = cj[3];
                  cj_inf  = cj[4];
                }
                pitch      = geometry->GetMaxAngularCoord(iMarker, iSpan) - geometry->GetMinAngularCoord(iMarker,iSpan);
                theta      = geometry->turbovertex[iMarker][iSpan][iVertex]->GetRelAngularCoord();
                deltaTheta = geometry->turbovertex[iMarker][iSpan][iVertex]->GetDeltaAngularCoord();
                TwoPiThetaFreq_Pitch = 2*PI_NUMBER*freq*theta/pitch;

                expArg = complex<su2double>(cos(TwoPiThetaFreq_Pitch)) - I*complex<su2double>(sin(TwoPiThetaFreq_Pitch));
                if (freq != 0){
                  cktemp_out1 +=  cj_out1*expArg*deltaTheta/pitch;
                  cktemp_out2 +=  cj_out2*expArg*deltaTheta/pitch;
                  cktemp_inf  +=  cj_inf*expArg*deltaTheta/pitch;
                }
                else{
                  cktemp_inf += complex<su2double>(0.0,0.0);
                  cktemp_out1 += complex<su2double>(0.0,0.0);
                  cktemp_out2 += complex<su2double>(0.0,0.0);
                }
              }

            }
          }
        }
      }

#ifdef HAVE_MPI
      MyRe_inf = cktemp_inf.real(); Re_inf = 0.0;
      MyIm_inf = cktemp_inf.imag(); Im_inf = 0.0;
      cktemp_inf = complex<su2double>(0.0,0.0);

      MyRe_out1 = cktemp_out1.real(); Re_out1 = 0.0;
      MyIm_out1 = cktemp_out1.imag(); Im_out1 = 0.0;
      cktemp_out1 = complex<su2double>(0.0,0.0);

      MyRe_out2 = cktemp_out2.real(); Re_out2 = 0.0;
      MyIm_out2 = cktemp_out2.imag(); Im_out2 = 0.0;
      cktemp_out2 = complex<su2double>(0.0,0.0);


      SU2_MPI::Allreduce(&MyRe_inf, &Re_inf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&MyIm_inf, &Im_inf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&MyRe_out1, &Re_out1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&MyIm_out1, &Im_out1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&MyRe_out2, &Re_out2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&MyIm_out2, &Im_out2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      cktemp_inf = complex<su2double>(Re_inf,Im_inf);
      cktemp_out1 = complex<su2double>(Re_out1,Im_out1);
      cktemp_out2 = complex<su2double>(Re_out2,Im_out2);

#endif

      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
          if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
            if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
              /*-----this is only valid 2D ----*/
              if (marker_flag == INFLOW){
                CkInflow[iMarker][iSpan][k]= cktemp_inf;
              }else{
                CkOutflow1[iMarker][iSpan][k]=cktemp_out1;
                CkOutflow2[iMarker][iSpan][k]=cktemp_out2;
              }
            }
          }
        }
      }
    }
  }

  delete [] turboVelocity;
  delete [] turboNormal;
  delete [] Velocity_i;
  delete [] deltaprim;
  delete [] cj;

}

void CEulerSolver::BC_Giles(CGeometry *geometry, CSolver **solver_container,
    CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iVar, jVar, iSpan;
  unsigned long  iPoint, Point_Normal, oldVertex, k, kend, kend_max;
  long iVertex;
  su2double  *UnitNormal, *turboVelocity, *turboNormal;

  su2double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, Density_b, Pressure_b, Temperature_b;
  su2double *Velocity_i, Velocity2_i, Energy_i, StaticEnergy_i, Density_i, Pressure_i;
  su2double Pressure_e;
  su2double *V_boundary, *V_domain, *S_boundary, *S_domain;
  unsigned short  iZone     = config->GetiZone();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  unsigned short nSpanWiseSections = geometry->GetnSpanWiseSections(config->GetMarker_All_TurbomachineryFlag(val_marker));
  su2double relfacAvgCfg       = config->GetGiles_RelaxFactorAverage(Marker_Tag);
  su2double relfacFouCfg       = config->GetGiles_RelaxFactorFourier(Marker_Tag);
  su2double *Normal;
  su2double TwoPiThetaFreq_Pitch, pitch,theta;
  su2double *SpanWiseValues = NULL;
  su2double spanPercent, extrarelfacAvg = 0.0, deltaSpan = 0.0, relfacAvg, relfacFou, coeffrelfacAvg = 0.0;
  unsigned short Turbo_Flag;

  Normal                = new su2double[nDim];
  turboNormal 	        = new su2double[nDim];
  UnitNormal            = new su2double[nDim];
  turboVelocity         = new su2double[nDim];
  Velocity_i            = new su2double[nDim];
  Velocity_b            = new su2double[nDim];


  su2double AverageSoundSpeed, *AverageTurboMach, AverageEntropy, AverageEnthalpy;
  AverageTurboMach = new su2double[nDim];
  S_boundary       = new su2double[8];

  su2double  AvgMach , *cj, GilesBeta, *delta_c, **R_Matrix, *deltaprim, **R_c_inv,**R_c, alphaIn_BC, gammaIn_BC = 0,
      P_Total, T_Total, *FlowDir, Enthalpy_BC, Entropy_BC, *R, *c_avg,*dcjs, Beta_inf2, c2js_Re, c3js_Re, cOutjs_Re, avgVel2 =0.0;

  long freq;

  delta_c       = new su2double[nVar];
  deltaprim     = new su2double[nVar];
  cj            = new su2double[nVar];
  R_Matrix      = new su2double*[nVar];
  R_c           = new su2double*[nVar-1];
  R_c_inv       = new su2double*[nVar-1];
  R             = new su2double[nVar-1];
  c_avg         = new su2double[nVar];
  dcjs          = new su2double[nVar];

  for (iVar = 0; iVar < nVar; iVar++)
  {
    R_Matrix[iVar] = new su2double[nVar];
    c_avg[iVar]    =  0.0;
    dcjs[iVar]     =  0.0;
  }
  for (iVar = 0; iVar < nVar-1; iVar++)
  {
    R_c[iVar] = new su2double[nVar-1];
    R_c_inv[iVar] = new su2double[nVar-1];
  }


  complex<su2double> I, c2ks, c2js, c3ks, c3js, c4ks, c4js, cOutks, cOutjs, Beta_inf;
  I = complex<su2double>(0.0,1.0);

  /*--- Compute coeff for under relaxation of Avg and Fourier Coefficient for hub and shroud---*/
  if (nDim == 3){
    extrarelfacAvg  = config->GetExtraRelFacGiles(0);
    spanPercent     = config->GetExtraRelFacGiles(1);
    Turbo_Flag      = config->GetMarker_All_TurbomachineryFlag(val_marker);
    SpanWiseValues  = geometry->GetSpanWiseValue(Turbo_Flag);
    deltaSpan       = SpanWiseValues[nSpanWiseSections-1]*spanPercent;
    coeffrelfacAvg  = (relfacAvgCfg - extrarelfacAvg)/deltaSpan;
  }

  for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){
    /*--- Compute under relaxation for the Hub and Shroud Avg and Fourier Coefficient---*/
    if(nDim == 3){
      if(SpanWiseValues[iSpan] <= SpanWiseValues[0] + deltaSpan){
        relfacAvg = extrarelfacAvg + coeffrelfacAvg*(SpanWiseValues[iSpan] - SpanWiseValues[0]);
        relfacFou = 0.0;
      }
      else if(SpanWiseValues[iSpan] >= SpanWiseValues[nSpanWiseSections -1] - deltaSpan){
        relfacAvg = extrarelfacAvg - coeffrelfacAvg*(SpanWiseValues[iSpan] - SpanWiseValues[nSpanWiseSections -1]);
        relfacFou = 0.0;
      }
      else{
        relfacAvg = relfacAvgCfg;
        relfacFou = relfacFouCfg;
      }
    }
    else{
      {
        relfacAvg = relfacAvgCfg;
        relfacFou = relfacFouCfg;
      }
    }

    FluidModel->SetTDState_Prho(AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan]);
    AverageSoundSpeed = FluidModel->GetSoundSpeed();
    AverageTurboMach[0] = AverageTurboVelocity[val_marker][iSpan][0]/AverageSoundSpeed;
    AverageTurboMach[1] = AverageTurboVelocity[val_marker][iSpan][1]/AverageSoundSpeed;

    if(grid_movement){
      AverageTurboMach[1] -= geometry->GetAverageTangGridVel(val_marker,iSpan)/AverageSoundSpeed;
    }

    AvgMach = AverageTurboMach[0]*AverageTurboMach[0] + AverageTurboMach[1]*AverageTurboMach[1];

    kend     = geometry->GetnFreqSpan(val_marker, iSpan);
    kend_max = geometry->GetnFreqSpanMax(config->GetMarker_All_TurbomachineryFlag(val_marker));
    conv_numerics->GetRMatrix(AverageSoundSpeed, AverageDensity[val_marker][iSpan], R_Matrix);

    switch(config->GetKind_Data_Giles(Marker_Tag)){

    case TOTAL_CONDITIONS_PT:

      /*--- Retrieve the specified total conditions for this inlet. ---*/
      P_Total  = config->GetGiles_Var1(Marker_Tag);
      T_Total  = config->GetGiles_Var2(Marker_Tag);
      FlowDir = config->GetGiles_FlowDir(Marker_Tag);
      alphaIn_BC = atan(FlowDir[1]/FlowDir[0]);

      gammaIn_BC = 0;
      if (nDim == 3){
        gammaIn_BC = FlowDir[2]; //atan(FlowDir[2]/FlowDir[0]);
      }

      /*--- Non-dim. the inputs---*/
      P_Total /= config->GetPressure_Ref();
      T_Total /= config->GetTemperature_Ref();

      /* --- Computes the total state --- */
      FluidModel->SetTDState_PT(P_Total, T_Total);
      Enthalpy_BC = FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
      Entropy_BC = FluidModel->GetEntropy();


      /* --- Computes the inverse matrix R_c --- */
      conv_numerics->ComputeResJacobianGiles(FluidModel, AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan], AverageTurboVelocity[val_marker][iSpan], alphaIn_BC, gammaIn_BC, R_c, R_c_inv);

      FluidModel->SetTDState_Prho(AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan]);
      AverageEnthalpy = FluidModel->GetStaticEnergy() + AveragePressure[val_marker][iSpan]/AverageDensity[val_marker][iSpan];
      AverageEntropy  = FluidModel->GetEntropy();

      avgVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) avgVel2 += AverageVelocity[val_marker][iSpan][iDim]*AverageVelocity[val_marker][iSpan][iDim];
      if (nDim == 2){
        R[0] = -(AverageEntropy - Entropy_BC);
        R[1] = -(AverageTurboVelocity[val_marker][iSpan][1] - tan(alphaIn_BC)*AverageTurboVelocity[val_marker][iSpan][0]);
        R[2] = -(AverageEnthalpy + 0.5*avgVel2 - Enthalpy_BC);
      }

      else{
        R[0] = -(AverageEntropy - Entropy_BC);
        R[1] = -(AverageTurboVelocity[val_marker][iSpan][1] - tan(alphaIn_BC)*AverageTurboVelocity[val_marker][iSpan][0]);
        R[2] = -(AverageTurboVelocity[val_marker][iSpan][2] - tan(gammaIn_BC)*AverageTurboVelocity[val_marker][iSpan][0]);
        R[3] = -(AverageEnthalpy + 0.5*avgVel2 - Enthalpy_BC);

      }
      /* --- Compute the avg component  c_avg = R_c^-1 * R --- */
      for (iVar = 0; iVar < nVar-1; iVar++){
        c_avg[iVar] = 0.0;
        for (jVar = 0; jVar < nVar-1; jVar++){
          c_avg[iVar] += R_c_inv[iVar][jVar]*R[jVar];
        }
      }
      break;

    case TOTAL_CONDITIONS_PT_1D:

      /*--- Retrieve the specified total conditions for this inlet. ---*/
      P_Total  = config->GetGiles_Var1(Marker_Tag);
      T_Total  = config->GetGiles_Var2(Marker_Tag);
      FlowDir = config->GetGiles_FlowDir(Marker_Tag);
      alphaIn_BC = atan(FlowDir[1]/FlowDir[0]);

      gammaIn_BC = 0;
      if (nDim == 3){
        // Review definition of angle
        gammaIn_BC = FlowDir[2]; //atan(FlowDir[2]/FlowDir[0]);
      }

      /*--- Non-dim. the inputs---*/
      P_Total /= config->GetPressure_Ref();
      T_Total /= config->GetTemperature_Ref();

      /* --- Computes the total state --- */
      FluidModel->SetTDState_PT(P_Total, T_Total);
      Enthalpy_BC = FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
      Entropy_BC = FluidModel->GetEntropy();


      /* --- Computes the inverse matrix R_c --- */
      conv_numerics->ComputeResJacobianGiles(FluidModel, AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan], AverageTurboVelocity[val_marker][iSpan], alphaIn_BC, gammaIn_BC, R_c, R_c_inv);

      FluidModel->SetTDState_Prho(AveragePressure[val_marker][nSpanWiseSections], AverageDensity[val_marker][nSpanWiseSections]);
      AverageEnthalpy = FluidModel->GetStaticEnergy() + AveragePressure[val_marker][nSpanWiseSections]/AverageDensity[val_marker][nSpanWiseSections];
      AverageEntropy  = FluidModel->GetEntropy();


      avgVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) avgVel2 += AverageVelocity[val_marker][iSpan][iDim]*AverageVelocity[val_marker][iSpan][iDim];
      if (nDim == 2){
        R[0] = -(AverageEntropy - Entropy_BC);
        R[1] = -(AverageTurboVelocity[val_marker][nSpanWiseSections][1] - tan(alphaIn_BC)*AverageTurboVelocity[val_marker][nSpanWiseSections][0]);
        R[2] = -(AverageEnthalpy + 0.5*avgVel2 - Enthalpy_BC);
      }

      else{
        R[0] = -(AverageEntropy - Entropy_BC);
        R[1] = -(AverageTurboVelocity[val_marker][nSpanWiseSections][1] - tan(alphaIn_BC)*AverageTurboVelocity[val_marker][nSpanWiseSections][0]);
        R[2] = -(AverageTurboVelocity[val_marker][nSpanWiseSections][2] - tan(gammaIn_BC)*AverageTurboVelocity[val_marker][nSpanWiseSections][0]);
        R[3] = -(AverageEnthalpy + 0.5*avgVel2 - Enthalpy_BC);

      }
      /* --- Compute the avg component  c_avg = R_c^-1 * R --- */
      for (iVar = 0; iVar < nVar-1; iVar++){
        c_avg[iVar] = 0.0;
        for (jVar = 0; jVar < nVar-1; jVar++){
          c_avg[iVar] += R_c_inv[iVar][jVar]*R[jVar];
        }
      }
      break;

    case MIXING_IN: case MIXING_OUT:

      /* --- Compute average jump of primitive at the mixing-plane interface--- */
      deltaprim[0] = ExtAverageDensity[val_marker][iSpan] - AverageDensity[val_marker][iSpan];
      deltaprim[1] = ExtAverageTurboVelocity[val_marker][iSpan][0] - AverageTurboVelocity[val_marker][iSpan][0];
      deltaprim[2] = ExtAverageTurboVelocity[val_marker][iSpan][1] - AverageTurboVelocity[val_marker][iSpan][1];
      if (nDim == 2){
        deltaprim[3] = ExtAveragePressure[val_marker][iSpan] - AveragePressure[val_marker][iSpan];
      }
      else
      {
        deltaprim[3] = ExtAverageTurboVelocity[val_marker][iSpan][2] - AverageTurboVelocity[val_marker][iSpan][2];
        deltaprim[4] = ExtAveragePressure[val_marker][iSpan] - AveragePressure[val_marker][iSpan];
      }


      /* --- Compute average jump of charachteristic variable at the mixing-plane interface--- */
      FluidModel->SetTDState_Prho(AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan]);
      AverageSoundSpeed = FluidModel->GetSoundSpeed();
      conv_numerics->GetCharJump(AverageSoundSpeed, AverageDensity[val_marker][iSpan], deltaprim, c_avg);
      break;

    case MIXING_IN_1D: case MIXING_OUT_1D:

      /* --- Compute average jump of primitive at the mixing-plane interface--- */
      deltaprim[0] = ExtAverageDensity[val_marker][nSpanWiseSections] - AverageDensity[val_marker][nSpanWiseSections];
      deltaprim[1] = ExtAverageTurboVelocity[val_marker][nSpanWiseSections][0] - AverageTurboVelocity[val_marker][nSpanWiseSections][0];
      deltaprim[2] = ExtAverageTurboVelocity[val_marker][nSpanWiseSections][1] - AverageTurboVelocity[val_marker][nSpanWiseSections][1];
      if (nDim == 2){
        deltaprim[3] = ExtAveragePressure[val_marker][nSpanWiseSections] - AveragePressure[val_marker][nSpanWiseSections];
      }
      else
      {
        deltaprim[3] = ExtAverageTurboVelocity[val_marker][nSpanWiseSections][2] - AverageTurboVelocity[val_marker][nSpanWiseSections][2];
        deltaprim[4] = ExtAveragePressure[val_marker][nSpanWiseSections] - AveragePressure[val_marker][nSpanWiseSections];
      }

      /* --- Compute average jump of charachteristic variable at the mixing-plane interface--- */
      FluidModel->SetTDState_Prho(AveragePressure[val_marker][nSpanWiseSections], AverageDensity[val_marker][nSpanWiseSections]);
      AverageSoundSpeed = FluidModel->GetSoundSpeed();
      conv_numerics->GetCharJump(AverageSoundSpeed, AverageDensity[val_marker][nSpanWiseSections], deltaprim, c_avg);
      break;


    case STATIC_PRESSURE:
      Pressure_e = config->GetGiles_Var1(Marker_Tag);
      Pressure_e /= config->GetPressure_Ref();

      /* --- Compute avg characteristic jump  --- */
      if (nDim == 2){
        c_avg[3] = -2.0*(AveragePressure[val_marker][iSpan]-Pressure_e);
      }
      else
      {
        c_avg[4] = -2.0*(AveragePressure[val_marker][iSpan]-Pressure_e);
      }
      break;

    case STATIC_PRESSURE_1D:
      Pressure_e = config->GetGiles_Var1(Marker_Tag);
      Pressure_e /= config->GetPressure_Ref();

      /* --- Compute avg characteristic jump  --- */
      if (nDim == 2){
        c_avg[3] = -2.0*(AveragePressure[val_marker][nSpanWiseSections]-Pressure_e);
      }
      else
      {
        c_avg[4] = -2.0*(AveragePressure[val_marker][nSpanWiseSections]-Pressure_e);
      }
      break;

    case RADIAL_EQUILIBRIUM:
      Pressure_e = RadialEquilibriumPressure[val_marker][iSpan];

      /* --- Compute avg characteristic jump  --- */
      c_avg[4] = -2.0*(AveragePressure[val_marker][iSpan]-Pressure_e);

      break;

    }

    /*--- Loop over all the vertices on this boundary marker ---*/

    for (iVertex = 0; iVertex < geometry->nVertexSpan[val_marker][iSpan]; iVertex++) {

      /*--- using the other vertex information for retrieving some information ---*/
      oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();
      V_boundary= GetCharacPrimVar(val_marker, oldVertex);

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention),
       * 		this normal is scaled with the area of the face of the element  ---*/
      geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- find the node related to the vertex ---*/
      iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- Normalize Normal vector for this vertex (already for outward convention) ---*/
      geometry->turbovertex[val_marker][iSpan][iVertex]->GetNormal(UnitNormal);
      geometry->turbovertex[val_marker][iSpan][iVertex]->GetTurboNormal(turboNormal);

      /*--- Retrieve solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();

      /*--- Retrieve domain Secondary variables ---*/
      S_domain = node[iPoint]->GetSecondary();


      /*--- Compute the internal state u_i ---*/
      Velocity2_i = 0;
      for (iDim = 0; iDim < nDim; iDim++)
      {
        Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
        Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
      }
      
      
      Density_i = node[iPoint]->GetDensity();
      
      Energy_i = node[iPoint]->GetEnergy();
      StaticEnergy_i = Energy_i - 0.5*Velocity2_i;
      
      FluidModel->SetTDState_rhoe(Density_i, StaticEnergy_i);
      
      Pressure_i = FluidModel->GetPressure();

      ComputeTurboVelocity(Velocity_i, turboNormal, turboVelocity, config->GetMarker_All_TurbomachineryFlag(val_marker),config->GetKind_TurboMachinery(iZone));
      if (nDim == 2){
        deltaprim[0] = Density_i - AverageDensity[val_marker][iSpan];
        deltaprim[1] = turboVelocity[0] - AverageTurboVelocity[val_marker][iSpan][0];
        deltaprim[2] = turboVelocity[1] - AverageTurboVelocity[val_marker][iSpan][1];
        deltaprim[3] = Pressure_i - AveragePressure[val_marker][iSpan];
      }
      else{
        deltaprim[0] = Density_i - AverageDensity[val_marker][iSpan];
        deltaprim[1] = turboVelocity[0] - AverageTurboVelocity[val_marker][iSpan][0];
        deltaprim[2] = turboVelocity[1] - AverageTurboVelocity[val_marker][iSpan][1];
        deltaprim[3] = turboVelocity[2] - AverageTurboVelocity[val_marker][iSpan][2];
        deltaprim[4] = Pressure_i - AveragePressure[val_marker][iSpan];
      }

      FluidModel->SetTDState_Prho(AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan]);
      AverageSoundSpeed = FluidModel->GetSoundSpeed();
      conv_numerics->GetCharJump(AverageSoundSpeed, AverageDensity[val_marker][iSpan], deltaprim, cj);


      pitch      = geometry->GetMaxAngularCoord(val_marker, iSpan) - geometry->GetMinAngularCoord(val_marker,iSpan);
      theta      = geometry->turbovertex[val_marker][iSpan][iVertex]->GetRelAngularCoord();

      switch(config->GetKind_Data_Giles(Marker_Tag))
      {

      //Done, generilize for 3D case
      //TODO(turbo), generilize for Inlet and Outlet in for backflow treatment

      case TOTAL_CONDITIONS_PT: case MIXING_IN:case TOTAL_CONDITIONS_PT_1D: case MIXING_IN_1D:
        if(config->GetSpatialFourier()){
          if (AvgMach <= 1.0){
            Beta_inf= I*complex<su2double>(sqrt(1.0 - AvgMach));
            c2js = complex<su2double>(0.0,0.0);
            c3js = complex<su2double>(0.0,0.0);
            for(k=0; k < 2*kend_max+1; k++){
              freq = k - kend_max;
              if(freq >= (long)(-kend) && freq <= (long)(kend) && AverageTurboMach[0] > config->GetAverageMachLimit()){
                TwoPiThetaFreq_Pitch = 2*PI_NUMBER*freq*theta/pitch;

                c2ks = -CkInflow[val_marker][iSpan][k]*complex<su2double>(Beta_inf + AverageTurboMach[1])/complex<su2double>( 1.0 + AverageTurboMach[0]);
                c3ks =  CkInflow[val_marker][iSpan][k]*complex<su2double>(Beta_inf + AverageTurboMach[1])/complex<su2double>( 1.0 + AverageTurboMach[0]);
                c3ks *= complex<su2double>(Beta_inf + AverageTurboMach[1])/complex<su2double>( 1.0 + AverageTurboMach[0]);
                c2js += c2ks*(complex<su2double>(cos(TwoPiThetaFreq_Pitch))+I*complex<su2double>(sin(TwoPiThetaFreq_Pitch)));
                c3js += c3ks*(complex<su2double>(cos(TwoPiThetaFreq_Pitch))+I*complex<su2double>(sin(TwoPiThetaFreq_Pitch)));
              }
              else{
                c2js += complex<su2double>(0.0,0.0);
                c3js += complex<su2double>(0.0,0.0);
              }
            }
            c2js_Re = c2js.real();
            c3js_Re = c3js.real();

            if (nDim == 2){
              dcjs[0] = 0.0     - cj[0];
              dcjs[1] = c2js_Re - cj[1];
              dcjs[2] = c3js_Re - cj[2];
            }else{
              dcjs[0] = 0.0     - cj[0];
              dcjs[1] = c2js_Re - cj[1];
              dcjs[2] = 0.0     - cj[2];
              dcjs[3] = c3js_Re - cj[3];
            }

          }else{
            if (AverageTurboVelocity[val_marker][iSpan][1] >= 0.0){
              Beta_inf2= -sqrt(AvgMach - 1.0);
            }else{
              Beta_inf2= sqrt(AvgMach-1.0);
            }
            if (nDim == 2){
              c2js_Re = -cj[3]*(Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
              c3js_Re = cj[3]*(Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
              c3js_Re *= (Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
            }else{
              c2js_Re = -cj[4]*(Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
              c3js_Re = cj[4]*(Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
              c3js_Re *= (Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
            }


            if (nDim == 2){
              dcjs[0] = 0.0     - cj[0];
              dcjs[1] = c2js_Re - cj[1];
              dcjs[2] = c3js_Re - cj[2];
            }else{
              dcjs[0] = 0.0     - cj[0];
              dcjs[1] = c2js_Re - cj[1];
              dcjs[2] = 0.0     - cj[2];
              dcjs[3] = c3js_Re - cj[3];
            }
          }
        }
        else{
          if (nDim == 2){
            dcjs[0] = 0.0;
            dcjs[1] = 0.0;
            dcjs[2] = 0.0;
          }else{
            dcjs[0] = 0.0;
            dcjs[1] = 0.0;
            dcjs[2] = 0.0;
            dcjs[3] = 0.0;
          }
        }
        /* --- Impose Inlet BC Reflecting--- */
        delta_c[0] = relfacAvg*c_avg[0] + relfacFou*dcjs[0];
        delta_c[1] = relfacAvg*c_avg[1] + relfacFou*dcjs[1];
        delta_c[2] = relfacAvg*c_avg[2] + relfacFou*dcjs[2];
        if (nDim == 2){
          delta_c[3] = cj[3];
        }else{
          delta_c[3] = relfacAvg*c_avg[3] + relfacFou*dcjs[3];
          delta_c[4] = cj[4];
        }
        break;


      case STATIC_PRESSURE:case STATIC_PRESSURE_1D:case MIXING_OUT:case RADIAL_EQUILIBRIUM:case MIXING_OUT_1D:

        /* --- implementation of Giles BC---*/
        if(config->GetSpatialFourier()){
          if (AvgMach > 1.0){
            /* --- supersonic Giles implementation ---*/
            if (AverageTurboVelocity[val_marker][iSpan][1] >= 0.0){
              GilesBeta= -sqrt(AvgMach - 1.0);

            }else{
              GilesBeta= sqrt(AvgMach - 1.0);
            }
            if(nDim == 2){
              cOutjs_Re= (2.0 * AverageTurboMach[0])/(GilesBeta - AverageTurboMach[1])*cj[1] - (GilesBeta + AverageTurboMach[1])/(GilesBeta - AverageTurboMach[1])*cj[2];
            }
            else{
              cOutjs_Re= (2.0 * AverageTurboMach[0])/(GilesBeta - AverageTurboMach[1])*cj[1] - (GilesBeta + AverageTurboMach[1])/(GilesBeta - AverageTurboMach[1])*cj[3];
            }
            if (nDim == 2){
              dcjs[3] = cOutjs_Re - cj[3];
            }
            else{
              dcjs[4] = cOutjs_Re - cj[4];
            }
          }else{

            /* --- subsonic Giles implementation ---*/
            Beta_inf= I*complex<su2double>(sqrt(1.0  - AvgMach));
            cOutjs 	= complex<su2double>(0.0,0.0);
            for(k=0; k < 2*kend_max+1; k++){
              freq = k - kend_max;
              if(freq >= (long)(-kend) && freq <= (long)(kend) && AverageTurboMach[0] > config->GetAverageMachLimit()){
                TwoPiThetaFreq_Pitch = 2*PI_NUMBER*freq*theta/pitch;
                cOutks  = complex<su2double>(2.0 * AverageTurboMach[0])/complex<su2double>(Beta_inf - AverageTurboMach[1])*CkOutflow1[val_marker][iSpan][k];
                cOutks -= complex<su2double>(Beta_inf + AverageTurboMach[1])/complex<su2double>(Beta_inf - AverageTurboMach[1])*CkOutflow2[val_marker][iSpan][k];

                cOutjs += cOutks*(complex<su2double>(cos(TwoPiThetaFreq_Pitch)) + I*complex<su2double>(sin(TwoPiThetaFreq_Pitch)));
              }
              else{
                cOutjs +=complex<su2double>(0.0,0.0);
              }
            }
            cOutjs_Re = cOutjs.real();

            if (nDim == 2){
              dcjs[3] = cOutjs_Re - cj[3];
            }
            else{
              dcjs[4] = cOutjs_Re - cj[4];
            }
          }
        }
        else{
          if (nDim == 2){
            dcjs[3] = 0.0;
          }
          else{
            dcjs[4] = 0.0;
          }
        }
        /* --- Impose Outlet BC Non-Reflecting  --- */
        delta_c[0] = cj[0];
        delta_c[1] = cj[1];
        delta_c[2] = cj[2];
        if (nDim == 2){
          delta_c[3] = relfacAvg*c_avg[3] + relfacFou*dcjs[3];
        }
        else{
          delta_c[3] = cj[3];
          delta_c[4] = relfacAvg*c_avg[4] + relfacFou*dcjs[4];
        }


        /*--- Automatically impose supersonic autoflow ---*/
        if (abs(AverageTurboMach[0]) > 1.0000){
          delta_c[0] = 0.0;
          delta_c[1] = 0.0;
          delta_c[2] = 0.0;
          delta_c[2] = 0.0;
          if (nDim == 3)delta_c[4] = 0.0;
        }

        break;

      default:
        SU2_MPI::Error("Invalid Giles input!", CURRENT_FUNCTION);
        break;
      }
      
      /*--- Compute primitive jump from characteristic variables  ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        deltaprim[iVar]=0.0;
        for (jVar = 0; jVar < nVar; jVar++)
        {
          deltaprim[iVar] +=  R_Matrix[iVar][jVar]*delta_c[jVar];
        }
      }

      /*--- retrieve boundary variables ---*/
      Density_b = AverageDensity[val_marker][iSpan] + deltaprim[0];
      turboVelocity[0] = AverageTurboVelocity[val_marker][iSpan][0] + deltaprim[1];
      turboVelocity[1] = AverageTurboVelocity[val_marker][iSpan][1] + deltaprim[2];
      if(nDim == 2){
        Pressure_b = AveragePressure[val_marker][iSpan] + deltaprim[3];
      }
      else{
        turboVelocity[2] = AverageTurboVelocity[val_marker][iSpan][2] + deltaprim[3];
        Pressure_b = AveragePressure[val_marker][iSpan] + deltaprim[4];
      }


      ComputeBackVelocity(turboVelocity, turboNormal, Velocity_b, config->GetMarker_All_TurbomachineryFlag(val_marker), config->GetKind_TurboMachinery(iZone));
      Velocity2_b = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity2_b+= Velocity_b[iDim]*Velocity_b[iDim];
      }

      if(Pressure_b <= 0.0 || Density_b <= 0.0 ){
        Pressure_b = Pressure_i;
        Density_b = Density_i;
        Velocity2_b = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity_b[iDim] = Velocity_i[iDim];
          Velocity2_b+= Velocity_b[iDim]*Velocity_b[iDim];
        }
      }

      FluidModel->SetTDState_Prho(Pressure_b, Density_b);
      Energy_b = FluidModel->GetStaticEnergy() + 0.5*Velocity2_b;
      Temperature_b= FluidModel->GetTemperature();
      Enthalpy_b = Energy_b + Pressure_b/Density_b;

      /*--- Primitive variables, using the derived quantities ---*/
      V_boundary[0] = Temperature_b;
      for (iDim = 0; iDim < nDim; iDim++)
        V_boundary[iDim+1] = Velocity_b[iDim];
      V_boundary[nDim+1] = Pressure_b;
      V_boundary[nDim+2] = Density_b;
      V_boundary[nDim+3] = Enthalpy_b;

      S_boundary[0]= FluidModel->GetdPdrho_e();
      S_boundary[1]= FluidModel->GetdPde_rho();



      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_boundary);
      conv_numerics->SetSecondary(S_domain, S_boundary);


      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      if (viscous) {

        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_boundary[nDim+5] = FluidModel->GetLaminarViscosity();
        V_boundary[nDim+6] = node[iPoint]->GetEddyViscosity();
        V_boundary[nDim+7] = FluidModel->GetThermalConductivity();
        V_boundary[nDim+8] = FluidModel->GetCp();
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_boundary);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());


        /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/
        
        S_boundary[0]= FluidModel->GetdPdrho_e();
        S_boundary[1]= FluidModel->GetdPde_rho();
        
        S_boundary[2]= FluidModel->GetdTdrho_e();
        S_boundary[3]= FluidModel->GetdTde_rho();
        
        /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/
        
        S_boundary[4]= FluidModel->Getdmudrho_T();
        S_boundary[5]= FluidModel->GetdmudT_rho();
        
        S_boundary[6]= FluidModel->Getdktdrho_T();
        S_boundary[7]= FluidModel->GetdktdT_rho();
        
        visc_numerics->SetSecondary(S_domain, S_boundary);
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
  delete [] Velocity_b;
  delete [] Velocity_i;
  
  delete [] S_boundary;
  delete []	delta_c;
  delete []	deltaprim;
  delete []	cj;
  for (iVar = 0; iVar < nVar; iVar++)
  {
    delete [] R_Matrix[iVar];
  }
  for (iVar = 0; iVar < nVar-1; iVar++)
  {
    delete [] R_c_inv[iVar];
    delete [] R_c[iVar];
  }
  delete [] R_Matrix;
  delete [] R_c;
  delete [] R_c_inv;
  delete [] R;
  delete [] c_avg;
  delete [] dcjs;

  delete [] AverageTurboMach;
  delete [] UnitNormal;
  delete [] turboNormal;
  delete [] turboVelocity;
}

void CEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
  Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
  alpha, aa, bb, cc, dd, Area, UnitNormal[3];
  su2double *V_inlet, *V_domain;
  
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  su2double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
  su2double Gas_Constant       = config->GetGas_ConstantND();
  unsigned short Kind_Inlet = config->GetKind_Inlet();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];
    
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the inlet ---*/
    
    V_inlet = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Build the fictitious intlet state based on characteristics ---*/
      

      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
         therefore we can specify all but one state variable at the inlet.
         The outgoing Riemann invariant provides the final piece of info.
         Adapted from an original implementation in the Stanford University
         multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
         written by Edwin van der Weide, last modified 04-20-2009. ---*/

      switch (Kind_Inlet) {

        /*--- Total properties have been specified at the inlet. ---*/

        case TOTAL_CONDITIONS:

          /*--- Retrieve the specified total conditions for this inlet. ---*/

          if (gravity) P_Total = Inlet_Ptotal[val_marker][iVertex] - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;
          else P_Total  = Inlet_Ptotal[val_marker][iVertex];
          T_Total  = Inlet_Ttotal[val_marker][iVertex];
          Flow_Dir = Inlet_FlowDir[val_marker][iVertex];

          /*--- Non-dim. the inputs if necessary. ---*/

          P_Total /= config->GetPressure_Ref();
          T_Total /= config->GetTemperature_Ref();

          /*--- Store primitives and set some variables for clarity. ---*/

          Density = V_domain[nDim+2];
          Velocity2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim] = V_domain[iDim+1];
            Velocity2 += Velocity[iDim]*Velocity[iDim];
          }
          Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
          Pressure    = V_domain[nDim+1];
          H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
          SoundSpeed2 = Gamma*Pressure/Density;

          /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/

          Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
          for (iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];

          /*--- Total speed of sound ---*/

          SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

          /*--- Dot product of normal and flow direction. This should
             be negative due to outward facing boundary normal convention. ---*/

          alpha = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            alpha += UnitNormal[iDim]*Flow_Dir[iDim];

          /*--- Coefficients in the quadratic equation for the velocity ---*/

          aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
          bb = -1.0*Gamma_Minus_One*alpha*Riemann;
          cc =  0.5*Gamma_Minus_One*Riemann*Riemann
              -2.0*SoundSpeed_Total2/Gamma_Minus_One;

          /*--- Solve quadratic equation for velocity magnitude. Value must
             be positive, so the choice of root is clear. ---*/

          dd = bb*bb - 4.0*aa*cc;
          dd = sqrt(max(0.0, dd));
          Vel_Mag   = (-bb + dd)/(2.0*aa);
          Vel_Mag   = max(0.0, Vel_Mag);
          Velocity2 = Vel_Mag*Vel_Mag;

          /*--- Compute speed of sound from total speed of sound eqn. ---*/

          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

          /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/

          Mach2 = Velocity2/SoundSpeed2;
          Mach2 = min(1.0, Mach2);
          Velocity2   = Mach2*SoundSpeed2;
          Vel_Mag     = sqrt(Velocity2);
          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

          /*--- Compute new velocity vector at the inlet ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

          /*--- Static temperature from the speed of sound relation ---*/

          Temperature = SoundSpeed2/(Gamma*Gas_Constant);

          /*--- Static pressure using isentropic relation at a point ---*/

          Pressure = P_Total*pow((Temperature/T_Total), Gamma/Gamma_Minus_One);

          /*--- Density at the inlet from the gas law ---*/

          Density = Pressure/(Gas_Constant*Temperature);

          /*--- Using pressure, density, & velocity, compute the energy ---*/

          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
          if (tkeNeeded) Energy += GetTke_Inf();

          /*--- Primitive variables, using the derived quantities ---*/

          V_inlet[0] = Temperature;
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Velocity[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;

          break;

          /*--- Mass flow has been specified at the inlet. ---*/

        case MASS_FLOW:

          /*--- Retrieve the specified mass flow for the inlet. ---*/

          Density  = Inlet_Ttotal[val_marker][iVertex];
          Vel_Mag  = Inlet_Ptotal[val_marker][iVertex];
          Flow_Dir = Inlet_FlowDir[val_marker][iVertex];

          /*--- Non-dim. the inputs if necessary. ---*/

          Density /= config->GetDensity_Ref();
          Vel_Mag /= config->GetVelocity_Ref();

          /*--- Get primitives from current inlet state. ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = node[iPoint]->GetVelocity(iDim);
          Pressure    = node[iPoint]->GetPressure();
          SoundSpeed2 = Gamma*Pressure/V_domain[nDim+2];

          /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/

          Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
          for (iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];

          /*--- Speed of sound squared for fictitious inlet state ---*/

          SoundSpeed2 = Riemann;
          for (iDim = 0; iDim < nDim; iDim++)
            SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitNormal[iDim];

          SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
          SoundSpeed2 = SoundSpeed2*SoundSpeed2;

          /*--- Pressure for the fictitious inlet state ---*/

          Pressure = SoundSpeed2*Density/Gamma;

          /*--- Energy for the fictitious inlet state ---*/

          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Vel_Mag*Vel_Mag;
          if (tkeNeeded) Energy += GetTke_Inf();

          /*--- Primitive variables, using the derived quantities ---*/

          V_inlet[0] = Pressure / ( Gas_Constant * Density);
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;

          break;
      }
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      if (viscous) {
//
//        /*--- Set laminar and eddy viscosity at the infinity ---*/
//
//        V_inlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
//        V_inlet[nDim+6] = node[iPoint]->GetEddyViscosity();
//
//        /*--- Set the normal vector and the coordinates ---*/
//
//        visc_numerics->SetNormal(Normal);
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//
//        /*--- Primitive variables, and gradient ---*/
//
//        visc_numerics->SetPrimitive(V_domain, V_inlet);
//        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
//
//        /*--- Turbulent kinetic energy ---*/
//
//        if (config->GetKind_Turb_Model() == SST)
//          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
//
//        /*--- Compute and update residual ---*/
//
//        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.SubtractBlock(iPoint, Residual);
//
//        /*--- Jacobian contribution for implicit integration ---*/
//
//        if (implicit)
//          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//
//      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint;
  su2double Pressure, P_Exit, Velocity[3],
  Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
  Area, UnitNormal[3];
  su2double *V_outlet, *V_domain;
  
  bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  su2double Gas_Constant     = config->GetGas_ConstantND();
  bool grid_movement      = config->GetGrid_Movement();
  string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the outlet ---*/
    V_outlet = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Build the fictitious intlet state based on characteristics ---*/

      /*--- Retrieve the specified back pressure for this outlet. ---*/
      if (gravity) P_Exit = config->GetOutlet_Pressure(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;
      else P_Exit = config->GetOutlet_Pressure(Marker_Tag);

      /*--- Non-dim. the inputs if necessary. ---*/
      P_Exit = P_Exit/config->GetPressure_Ref();

      /*--- Check whether the flow is supersonic at the exit. The type
         of boundary update depends on this. ---*/
      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      Pressure   = V_domain[nDim+1];
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Mach_Exit  = sqrt(Velocity2)/SoundSpeed;

      if (Mach_Exit >= 1.0) {

        /*--- Supersonic exit flow: there are no incoming characteristics,
           so no boundary condition is necessary. Set outlet state to current
           state so that upwinding handles the direction of propagation. ---*/
        for (iVar = 0; iVar < nPrimVar; iVar++) V_outlet[iVar] = V_domain[iVar];

      } else {

        /*--- Subsonic exit flow: there is one incoming characteristic,
           therefore one variable can be specified (back pressure) and is used
           to update the conservative variables. Compute the entropy and the
           acoustic Riemann variable. These invariants, as well as the
           tangential velocity components, are extrapolated. Adapted from an
           original implementation in the Stanford University multi-block
           (SUmb) solver in the routine bcSubsonicOutflow.f90 by Edwin van
           der Weide, last modified 09-10-2007. ---*/

        Entropy = Pressure*pow(1.0/Density, Gamma);
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

        /*--- Compute the new fictious state at the outlet ---*/
        Density    = pow(P_Exit/Entropy,1.0/Gamma);
        Pressure   = P_Exit;
        SoundSpeed = sqrt(Gamma*P_Exit/Density);
        Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        Velocity2  = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }
        Energy = P_Exit/(Density*Gamma_Minus_One) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();

        /*--- Conservative variables, using the derived quantities ---*/
        V_outlet[0] = Pressure / ( Gas_Constant * Density);
        for (iDim = 0; iDim < nDim; iDim++)
          V_outlet[iDim+1] = Velocity[iDim];
        V_outlet[nDim+1] = Pressure;
        V_outlet[nDim+2] = Density;
        V_outlet[nDim+3] = Energy + Pressure/Density;

      }
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
//      /*--- Viscous contribution, commented out because serious convergence problems  ---*/
//
//      if (viscous) {
//
//        /*--- Set laminar and eddy viscosity at the infinity ---*/
//        V_outlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
//        V_outlet[nDim+6] = node[iPoint]->GetEddyViscosity();
//
//        /*--- Set the normal vector and the coordinates ---*/
//        visc_numerics->SetNormal(Normal);
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//
//        /*--- Primitive variables, and gradient ---*/
//        visc_numerics->SetPrimitive(V_domain, V_outlet);
//        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
//
//        /*--- Turbulent kinetic energy ---*/
//        if (config->GetKind_Turb_Model() == SST)
//          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
//
//        /*--- Compute and update residual ---*/
//        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.SubtractBlock(iPoint, Residual);
//
//        /*--- Jacobian contribution for implicit integration ---*/
//        if (implicit)
//         Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//
//      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
}

void CEulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                                       CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_inlet, *V_domain;
  
  su2double Density, Pressure, Temperature, Energy, *Velocity, Velocity2;
  su2double Gas_Constant = config->GetGas_ConstantND();
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];
  
  /*--- Supersonic inlet flow: there are no outgoing characteristics,
   so all flow variables can be imposed at the inlet.
   First, retrieve the specified values for the primitive variables. ---*/
  
  Temperature = config->GetInlet_Temperature(Marker_Tag);
  Pressure    = config->GetInlet_Pressure(Marker_Tag);
  Velocity    = config->GetInlet_Velocity(Marker_Tag);
  
  /*--- Density at the inlet from the gas law ---*/
  
  Density = Pressure/(Gas_Constant*Temperature);
  
  /*--- Non-dim. the inputs if necessary. ---*/
  
  Temperature = Temperature/config->GetTemperature_Ref();
  Pressure    = Pressure/config->GetPressure_Ref();
  Density     = Density/config->GetDensity_Ref();
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity[iDim] = Velocity[iDim]/config->GetVelocity_Ref();
  
  /*--- Compute the energy from the specified state ---*/
  
  Velocity2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity2 += Velocity[iDim]*Velocity[iDim];
  Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;
  if (tkeNeeded) Energy += GetTke_Inf();
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the outlet ---*/
    
    V_inlet = GetCharacPrimVar(val_marker, iVertex);
    
    /*--- Primitive variables, using the derived quantities ---*/
    
    V_inlet[0] = Temperature;
    for (iDim = 0; iDim < nDim; iDim++)
      V_inlet[iDim+1] = Velocity[iDim];
    V_inlet[nDim+1] = Pressure;
    V_inlet[nDim+2] = Density;
    V_inlet[nDim+3] = Energy + Pressure/Density;
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      if (viscous) {
//
//        /*--- Set laminar and eddy viscosity at the infinity ---*/
//
//        V_inlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
//        V_inlet[nDim+6] = node[iPoint]->GetEddyViscosity();
//
//        /*--- Set the normal vector and the coordinates ---*/
//
//        visc_numerics->SetNormal(Normal);
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//
//        /*--- Primitive variables, and gradient ---*/
//
//        visc_numerics->SetPrimitive(V_domain, V_inlet);
//        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
//
//        /*--- Turbulent kinetic energy ---*/
//
//        if (config->GetKind_Turb_Model() == SST)
//          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
//
//        /*--- Compute and update residual ---*/
//
//        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.SubtractBlock(iPoint, Residual);
//
//        /*--- Jacobian contribution for implicit integration ---*/
//
//        if (implicit)
//          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container,
                                        CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_outlet, *V_domain;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  su2double *Normal = new su2double[nDim];
  
  /*--- Supersonic outlet flow: there are no ingoing characteristics,
   so all flow variables can should be interpolated from the domain. ---*/
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Allocate the value at the outlet ---*/
      
      V_outlet = GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Primitive variables, using the derived quantities ---*/
      
      V_outlet[0] = V_domain[0];
      for (iDim = 0; iDim < nDim; iDim++)
        V_outlet[iDim+1] = V_domain[iDim+1];
      V_outlet[nDim+1] = V_domain[nDim+1];
      V_outlet[nDim+2] = V_domain[nDim+2];
      V_outlet[nDim+3] = V_domain[nDim+3];
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
 //     /*--- Viscous contribution, commented out because serious convergence problems ---*/
 //
 //     if (viscous) {
 //
 //       /*--- Set laminar and eddy viscosity at the infinity ---*/
 //
 //       V_outlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
 //       V_outlet[nDim+6] = node[iPoint]->GetEddyViscosity();
 //
 //       /*--- Set the normal vector and the coordinates ---*/
 //
 //       visc_numerics->SetNormal(Normal);
 //       visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
 //
 //       /*--- Primitive variables, and gradient ---*/
 //
 //       visc_numerics->SetPrimitive(V_domain, V_outlet);
 //       visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
 //
 //       /*--- Turbulent kinetic energy ---*/
 //
 //       if (config->GetKind_Turb_Model() == SST)
 //         visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
//
//        /*--- Compute and update residual ---*/
//
//        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.SubtractBlock(iPoint, Residual);
//
//        /*--- Jacobian contribution for implicit integration ---*/
//
//        if (implicit)
//          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double Pressure, Inflow_Pressure = 0.0, Velocity[3], Velocity2, Entropy, Target_Inflow_MassFlow = 0.0, Target_Inflow_Mach = 0.0, Density, Energy,
  Riemann, Area, UnitNormal[3], Vn, SoundSpeed, Vn_Exit, Inflow_Pressure_inc, Inflow_Pressure_old, Inflow_Mach_old, Inflow_MassFlow_old;
  su2double *V_inflow, *V_domain;
  
  bool grid_movement  = config->GetGrid_Movement();

  su2double DampingFactor = config->GetDamp_Engine_Inflow();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  unsigned short Kind_Engine_Inflow = config->GetKind_Engine_Inflow();
  su2double Gas_Constant = config->GetGas_ConstantND();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double Baseline_Press = 0.75 * config->GetPressure_FreeStreamND();
  bool Engine_HalfModel = config->GetEngine_HalfModel();

  su2double *Normal = new su2double[nDim];
  
  
  if (Kind_Engine_Inflow == FAN_FACE_MACH) {
    
    /*--- Retrieve the specified target fan face mach at the nacelle. ---*/
    
    Target_Inflow_Mach = config->GetEngineInflow_Target(Marker_Tag);
    
    /*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/
    
    Inflow_Pressure_old = config->GetInflow_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
    Inflow_Mach_old = config->GetInflow_Mach(Marker_Tag);
    
    /*--- Compute the pressure increment (note that increasing pressure decreases flow speed) ---*/
    
    Inflow_Pressure_inc = - (1.0 - (Inflow_Mach_old/Target_Inflow_Mach)) * Baseline_Press;
    
    /*--- Estimate the new fan face pressure ---*/
    
    Inflow_Pressure = (1.0 - DampingFactor)*Inflow_Pressure_old + DampingFactor * (Inflow_Pressure_old + Inflow_Pressure_inc);
    
  }
  
  if (Kind_Engine_Inflow == FAN_FACE_MDOT) {
    
    /*--- Retrieve the specified target mass flow (non-dimensional) at the nacelle. ---*/
    
    Target_Inflow_MassFlow = config->GetEngineInflow_Target(Marker_Tag) / (config->GetDensity_Ref() * config->GetVelocity_Ref());
    
    if (config->GetSystemMeasurements() == US) Target_Inflow_MassFlow /= 32.174;
    
    if (Engine_HalfModel) Target_Inflow_MassFlow /= 2.0;

    /*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/
    
    Inflow_Pressure_old = config->GetInflow_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
    Inflow_MassFlow_old = config->GetInflow_MassFlow(Marker_Tag);  // same here... it is a non dimensional value
    
    /*--- Compute the pressure increment (note that increasing pressure decreases flow speed) ---*/
    
    Inflow_Pressure_inc = - (1.0 - (Inflow_MassFlow_old/Target_Inflow_MassFlow)) * Baseline_Press;
    
    /*--- Estimate the new fan face pressure ---*/
    
    Inflow_Pressure = (1.0 - DampingFactor)*Inflow_Pressure_old + DampingFactor * (Inflow_Pressure_old + Inflow_Pressure_inc);
    
  }
  
  /*--- No iterative scheme if we provide the static pressure ---*/
  
  if (Kind_Engine_Inflow == FAN_FACE_PRESSURE) {
    
    /*--- Retrieve the specified pressure (non-dimensional) at the nacelle. ---*/
    
    Inflow_Pressure = config->GetEngineInflow_Target(Marker_Tag) / config->GetPressure_Ref();
    
  }
  
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the outlet ---*/
    
    V_inflow = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Subsonic nacelle inflow: there is one incoming characteristic,
       therefore one variable can be specified (back pressure) and is used
       to update the conservative variables.
       
       Compute the entropy and the acoustic variable. These
       riemann invariants, as well as the tangential velocity components,
       are extrapolated. ---*/
      
      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      Pressure   = V_domain[nDim+1];
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Entropy = Pressure*pow(1.0/Density, Gamma);
      Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;
      
      /*--- Compute the new fictious state at the outlet ---*/
      
      Density    = pow(Inflow_Pressure/Entropy,1.0/Gamma);
      Pressure   = Inflow_Pressure;
      SoundSpeed = sqrt(Gamma*Inflow_Pressure/Density);
      Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
      Velocity2  = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      
      Energy = Inflow_Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
      if (tkeNeeded) Energy += GetTke_Inf();
      
      /*--- Conservative variables, using the derived quantities ---*/
      
      V_inflow[0] = Pressure / ( Gas_Constant * Density);
      for (iDim = 0; iDim < nDim; iDim++)
        V_inflow[iDim+1] = Velocity[iDim];
      V_inflow[nDim+1] = Pressure;
      V_inflow[nDim+2] = Density;
      V_inflow[nDim+3] = Energy + Pressure/Density;
      V_inflow[nDim+4] = SoundSpeed;
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_inflow);
      
      /*--- Set grid movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      if (viscous) {
//
//        /*--- Set laminar and eddy viscosity at the infinity ---*/
//
//        V_inflow[nDim+5] = node[iPoint]->GetLaminarViscosity();
//        V_inflow[nDim+6] = node[iPoint]->GetEddyViscosity();
//
//        /*--- Set the normal vector and the coordinates ---*/
//
//        visc_numerics->SetNormal(Normal);
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//
//        /*--- Primitive variables, and gradient ---*/
//
//        visc_numerics->SetPrimitive(V_domain, V_inflow);
//        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
//
//        /*--- Turbulent kinetic energy ---*/
//
//        if (config->GetKind_Turb_Model() == SST)
//          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
//
//        /*--- Compute and update residual ---*/
//
//        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.SubtractBlock(iPoint, Residual);
//
//        /*--- Jacobian contribution for implicit integration ---*/
//
//        if (implicit)
//          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//
//      }
      
    }
  }
  
  delete [] Normal;
  
}


void CEulerSolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double Exhaust_Pressure, Exhaust_Temperature, Velocity[3], Velocity2, H_Exhaust, Temperature, Riemann, Area, UnitNormal[3], Pressure, Density, Energy, Mach2, SoundSpeed2, SoundSpeed_Exhaust2, Vel_Mag, alpha, aa, bb, cc, dd, Flow_Dir[3];
  su2double *V_exhaust, *V_domain, Target_Exhaust_Pressure, Exhaust_Pressure_old, Exhaust_Pressure_inc;
  
  su2double Gas_Constant = config->GetGas_ConstantND();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double DampingFactor = config->GetDamp_Engine_Exhaust();
  su2double Baseline_Press = 0.75 * config->GetPressure_FreeStreamND();
  
  su2double *Normal = new su2double[nDim];
  
  /*--- Retrieve the specified exhaust pressure in the engine (non-dimensional). ---*/
  
  Target_Exhaust_Pressure = config->GetExhaust_Pressure_Target(Marker_Tag) / config->GetPressure_Ref();
  
  /*--- Retrieve the old exhaust pressure in the engine exhaust (this has been computed in a preprocessing). ---*/
  
  Exhaust_Pressure_old = config->GetExhaust_Pressure(Marker_Tag);
  
  /*--- Compute the Pressure increment ---*/
  
  Exhaust_Pressure_inc = (1.0 - (Exhaust_Pressure_old/Target_Exhaust_Pressure)) * Baseline_Press;
  
  /*--- Estimate the new exhaust pressure ---*/
  
  Exhaust_Pressure = (1.0 - DampingFactor) * Exhaust_Pressure_old + DampingFactor * (Exhaust_Pressure_old + Exhaust_Pressure_inc);
  
  /*--- The temperature is given (no iteration is required) ---*/
  
  Exhaust_Temperature  = config->GetExhaust_Temperature_Target(Marker_Tag);
  Exhaust_Temperature /= config->GetTemperature_Ref();
  
  /*--- The pressure is given (no iteration is required) ---*/
  
  Exhaust_Pressure  = config->GetExhaust_Pressure_Target(Marker_Tag);
  Exhaust_Pressure /= config->GetPressure_Ref();
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the exhaust ---*/
    
    V_exhaust = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
       therefore we can specify all but one state variable at the inlet.
       The outgoing Riemann invariant provides the final piece of info. ---*/
      
      /*--- Store primitives and set some variables for clarity. ---*/
      
      Density = V_domain[nDim+2];
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
      Pressure    = V_domain[nDim+1];
      H_Exhaust   = (Gamma*Gas_Constant/Gamma_Minus_One)*Exhaust_Temperature;
      SoundSpeed2 = Gamma*Pressure/Density;
      
      /*--- Compute the acoustic Riemann invariant that is extrapolated
       from the domain interior. ---*/
      
      Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
      for (iDim = 0; iDim < nDim; iDim++)
        Riemann += Velocity[iDim]*UnitNormal[iDim];
      
      /*--- Total speed of sound ---*/
      
      SoundSpeed_Exhaust2 = Gamma_Minus_One*(H_Exhaust - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
      
      /*--- The flow direction is defined by the surface normal ---*/
      
      for (iDim = 0; iDim < nDim; iDim++)
        Flow_Dir[iDim] = -UnitNormal[iDim];
      
      /*--- Dot product of normal and flow direction. This should
       be negative due to outward facing boundary normal convention. ---*/
      
      alpha = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        alpha += UnitNormal[iDim]*Flow_Dir[iDim];
      
      /*--- Coefficients in the quadratic equation for the velocity ---*/
      
      aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
      bb = -1.0*Gamma_Minus_One*alpha*Riemann;
      cc =  0.5*Gamma_Minus_One*Riemann*Riemann - 2.0*SoundSpeed_Exhaust2/Gamma_Minus_One;
      
      /*--- Solve quadratic equation for velocity magnitude. Value must
       be positive, so the choice of root is clear. ---*/
      
      dd      = bb*bb - 4.0*aa*cc;
      dd      = sqrt(max(0.0, dd));
      Vel_Mag = (-bb + dd)/(2.0*aa);
      
      if (Vel_Mag >= 0.0) {
        
        Velocity2 = Vel_Mag*Vel_Mag;
        
        /*--- Compute speed of sound from total speed of sound eqn. ---*/
        
        SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;
        Mach2       = Velocity2/SoundSpeed2;
        Velocity2   = Mach2*SoundSpeed2;
        Vel_Mag     = sqrt(Velocity2);
        SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;
        
        /*--- Compute new velocity vector at the inlet ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
        
        /*--- Static temperature from the speed of sound relation ---*/
        
        Temperature = SoundSpeed2/(Gamma*Gas_Constant);
        
        /*--- Static pressure using isentropic relation at a point ---*/
        
        Pressure = Exhaust_Pressure*pow((Temperature/Exhaust_Temperature), Gamma/Gamma_Minus_One);
        
        /*--- Density at the exhaust from the gas law ---*/
        
        Density = Pressure/(Gas_Constant*Temperature);
        
        /*--- Using pressure, density, & velocity, compute the energy ---*/
        
        Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();
        
        /*--- Primitive variables, using the derived quantities ---*/
        
        V_exhaust[0] = Temperature;
        for (iDim = 0; iDim < nDim; iDim++)
          V_exhaust[iDim+1] = Velocity[iDim];
        V_exhaust[nDim+1] = Pressure;
        V_exhaust[nDim+2] = Density;
        V_exhaust[nDim+3] = Energy + Pressure/Density;
        V_exhaust[nDim+4] = sqrt(SoundSpeed2);
        
      } 
      /*--- The flow goes in the wrong direction ---*/
      
      else {
        
        V_exhaust[0] = V_domain[0];
        for (iDim = 0; iDim < nDim; iDim++)
          V_exhaust[iDim+1] = V_domain[iDim+1];
        V_exhaust[nDim+1] = V_domain[nDim+1];
        V_exhaust[nDim+2] = V_domain[nDim+2];
        V_exhaust[nDim+3] = V_domain[nDim+3];
        V_exhaust[nDim+4] = V_domain[nDim+4];
        
      }
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_exhaust);
      
      /*--- Set grid movement ---*/

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      if (viscous) {
//
//        /*--- Set laminar and eddy viscosity at the infinity ---*/
//
//        V_exhaust[nDim+5] = node[iPoint]->GetLaminarViscosity();
//        V_exhaust[nDim+6] = node[iPoint]->GetEddyViscosity();
//
//        /*--- Set the normal vector and the coordinates ---*/
//
//        visc_numerics->SetNormal(Normal);
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//
//        /*--- Primitive variables, and gradient ---*/
//
//        visc_numerics->SetPrimitive(V_domain, V_exhaust);
//        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
//
//        /*--- Turbulent kinetic energy ---*/
//
//        if (config->GetKind_Turb_Model() == SST)
//          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
//
//        /*--- Compute and update residual ---*/
//
//        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.SubtractBlock(iPoint, Residual);
//
//        /*--- Jacobian contribution for implicit integration ---*/
//
//        if (implicit)
//          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//
//      }
      
    }
  }
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker) {
  
  /*--- Call the Euler residual ---*/
  
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);
  
}

void CEulerSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                         CConfig *config) {
  
  unsigned long iVertex, jVertex, iPoint, Point_Normal = 0;
  unsigned short iDim, iVar, iMarker, nDonorVertex;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool viscous       = config->GetViscous();
  
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  su2double *tmp_residual = new su2double[nVar];
  
  su2double weight;
  su2double P_static, rho_static;
   
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {

          nDonorVertex = GetnSlidingStates(iMarker, iVertex);
          
          /*--- Initialize Residual, this will serve to accumulate the average ---*/

          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = 0.0;

          /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/

          for (jVertex = 0; jVertex < nDonorVertex; jVertex++){

            Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
              PrimVar_j[iVar] = GetSlidingState(iMarker, iVertex, iVar, jVertex);
            }
            
            /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/

            weight = GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

            /*--- Set primitive variables ---*/

            conv_numerics->SetPrimitive( PrimVar_i, PrimVar_j );
          
            if( !( config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS ) ) {
              Secondary_i = node[iPoint]->GetSecondary();

              P_static   = PrimVar_j[nDim+1];
              rho_static = PrimVar_j[nDim+2];           
              FluidModel->SetTDState_Prho(P_static, rho_static);

              Secondary_j[0] = FluidModel->GetdPdrho_e();
              Secondary_j[1] = FluidModel->GetdPde_rho();  

              conv_numerics->SetSecondary(Secondary_i, Secondary_j);
            }

            /*--- Set the normal vector ---*/
 
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) 
              Normal[iDim] = -Normal[iDim];

            conv_numerics->SetNormal(Normal);

            if (grid_movement)
              conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the convective residual using an upwind scheme ---*/

            conv_numerics->ComputeResidual(tmp_residual, Jacobian_i, Jacobian_j, config);

            /*--- Accumulate the residuals to compute the average ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
              Residual[iVar] += weight*tmp_residual[iVar];

          }

          /*--- Add Residuals and Jacobians ---*/
  
          LinSysRes.AddBlock(iPoint, Residual);
          if (implicit) 
            Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

          if (viscous) {
            
            /*--- Initialize Residual, this will serve to accumulate the average ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
              Residual[iVar] = 0.0;
              
            /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/
            
            for (jVertex = 0; jVertex < nDonorVertex; jVertex++){
              PrimVar_j[nDim+5] = GetSlidingState(iMarker, iVertex, nDim+5, jVertex); 
              PrimVar_j[nDim+6] = GetSlidingState(iMarker, iVertex, nDim+6, jVertex); 

              /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/
              
              weight = GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);
              
              /*--- Set the normal vector and the coordinates ---*/

              visc_numerics->SetNormal(Normal);
              visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

              /*--- Primitive variables, and gradient ---*/

              visc_numerics->SetPrimitive(PrimVar_i, PrimVar_j);
              visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

              /*--- Turbulent kinetic energy ---*/

              if (config->GetKind_Turb_Model() == SST)
                visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));

              /*--- Compute and update residual ---*/

              visc_numerics->ComputeResidual(tmp_residual, Jacobian_i, Jacobian_j, config);
              
              /*--- Accumulate the residuals to compute the average ---*/
              
              for (iVar = 0; iVar < nVar; iVar++)
                Residual[iVar] += weight*tmp_residual[iVar];
            }
          
            LinSysRes.SubtractBlock(iPoint, Residual);
 
            /*--- Jacobian contribution for implicit integration ---*/

            if (implicit)
              Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
            
          }
        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] tmp_residual;
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
}

void CEulerSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, GlobalIndex_iPoint, GlobalIndex_jPoint;
  unsigned short iDim, iVar;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  /*--- Do the send process, by the moment we are sending each
   node individually, this must be changed ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
    GlobalIndex_jPoint = GetDonorGlobalIndex(val_marker, iVertex);
    
    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
      
      /*--- Store the solution for both points ---*/
      
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
        PrimVar_j[iVar] = GetDonorPrimVar(val_marker, iVertex, iVar);
      }
      
      /*--- Set Conservative Variables ---*/
      
      numerics->SetPrimitive(PrimVar_i, PrimVar_j);
      
      /*--- Set Normal ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      numerics->SetNormal(Normal);
      
      /*--- Compute the convective residual using an upwind scheme ---*/
      
      numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add Residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
    
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CEulerSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, GlobalIndex_iPoint, GlobalIndex_jPoint;
  unsigned short iDim, iVar;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  /*--- Do the send process, by the moment we are sending each
   node individually, this must be changed ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
    GlobalIndex_jPoint = GetDonorGlobalIndex(val_marker, iVertex);
    
    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
      
      /*--- Store the solution for both points ---*/
      
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
        PrimVar_j[iVar] = GetDonorPrimVar(val_marker, iVertex, iVar);
      }
      
      /*--- Set Conservative Variables ---*/
      
      numerics->SetPrimitive(PrimVar_i, PrimVar_j);
      
      /*--- Set Normal ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      numerics->SetNormal(Normal);
      
      /*--- Compute the convective residual using an upwind scheme ---*/
      
      numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add Residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
    
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CEulerSolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker) {
  
  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker, true);
  
}

void CEulerSolver::BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                     CConfig *config, unsigned short val_marker) {
  
  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker, false);
  
}

void CEulerSolver::BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                              CConfig *config, unsigned short val_marker, bool inlet_surface) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, GlobalIndex_donor, GlobalIndex;
  su2double Pressure, Velocity[3], Target_Press_Jump, Target_Temp_Jump,
  Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Vn_Inlet, Mach_Outlet,
  Area, UnitNormal[3], *V_outlet, *V_domain, *V_inlet, P_Total, T_Total, H_Total, Temperature,
  Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag, alpha, aa, bb, cc, dd;
  su2double Factor, P_static, T_static, SoS_outlet, Rho_outlet, Rho_inlet;
  su2double Vel_normal_inlet[3], Vel_tangent_inlet[3], Vel_inlet[3];
  su2double Vel_normal_outlet[3], Vel_tangent_outlet[3], Vel_outlet[3];
  su2double Vel_normal_inlet_, Vel_tangent_inlet_, Vel_inlet_;
  su2double Vel_normal_outlet_, Vel_outlet_;
  
  bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  su2double Gas_Constant  = config->GetGas_ConstantND();
  bool grid_movement      = config->GetGrid_Movement();
  bool tkeNeeded          = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                            (config->GetKind_Turb_Model() == SST));
  bool ratio              = (config->GetActDisk_Jump() == RATIO);
  su2double SecondaryFlow = config->GetSecondaryFlow_ActDisk();
  
  su2double *Normal = new su2double[nDim];
  su2double *Flow_Dir = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    if (inlet_surface) {
      V_inlet = GetCharacPrimVar(val_marker, iVertex);
      V_outlet = GetDonorPrimVar(val_marker, iVertex);
    }
    else {
      V_outlet = GetCharacPrimVar(val_marker, iVertex);
      V_inlet = GetDonorPrimVar(val_marker, iVertex);
    }
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
    GlobalIndex_donor = GetDonorGlobalIndex(val_marker, iVertex);
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if ((geometry->node[iPoint]->GetDomain()) &&
        (GlobalIndex != GlobalIndex_donor)) {
      
      /*--- Normal vector for this vertex (negative for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node and jumps values ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      Target_Press_Jump = GetActDisk_DeltaP(val_marker, iVertex);
      Target_Temp_Jump = GetActDisk_DeltaT(val_marker, iVertex);
      
      if (inlet_surface) {
        if (ratio) { P_static = V_outlet[nDim+1]/Target_Press_Jump; T_static = V_outlet[0]/Target_Temp_Jump; }
        else { P_static = V_outlet[nDim+1] - Target_Press_Jump; T_static = V_outlet[0] - Target_Temp_Jump; }
      }
      else {
        if (ratio) { P_static = V_inlet[nDim+1]*Target_Press_Jump; T_static = V_inlet[0]*Target_Temp_Jump; }
        else { P_static = V_inlet[nDim+1] + Target_Press_Jump; T_static = V_inlet[0] + Target_Temp_Jump; }
      }
      
      /*--- Subsonic inlet ---*/
      
      if (inlet_surface) {
        
        /*--- Build the fictitious intlet state based on characteristics.
         Retrieve the specified back pressure for this inlet ---*/
        
        Density = V_domain[nDim+2];
        Velocity2 = 0.0; Vn = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = V_domain[iDim+1];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
          Vn += Velocity[iDim]*UnitNormal[iDim];
        }
        Pressure   = V_domain[nDim+1];
        SoundSpeed = sqrt(Gamma*Pressure/Density);
        
        Entropy = Pressure*pow(1.0/Density, Gamma);
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;
        
        /*--- Compute the new fictious state at the outlet ---*/
        
        Pressure   = P_static;
        Density    = pow(Pressure/Entropy,1.0/Gamma);
        SoundSpeed = sqrt(Gamma*Pressure/Density);
        Vn_Inlet    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        
        Velocity2  = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = Velocity[iDim] + (Vn_Inlet-Vn)*UnitNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }
        Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();
        
        /*--- Conservative variables, using the derived quantities ---*/
        
        V_inlet[0] = Pressure / ( Gas_Constant * Density);
        for (iDim = 0; iDim < nDim; iDim++)
          V_inlet[iDim+1] = Velocity[iDim];
        V_inlet[nDim+1] = Pressure;
        V_inlet[nDim+2] = Density;
        V_inlet[nDim+3] = Energy + Pressure/Density;
        V_inlet[nDim+4] = SoundSpeed;
        conv_numerics->SetPrimitive(V_domain, V_inlet);
        
      }
      
      /*--- Subsonic outlet ---*/
      
      else {
        
        FluidModel->SetTDState_PT(P_static, T_static);
        SoS_outlet = FluidModel->GetSoundSpeed();
        Rho_outlet = FluidModel->GetDensity();
        
        /*--- We use the velocity and the density from the flow inlet
         to evaluate flow direction and mass flow ---*/
        
        Rho_inlet = V_inlet[nDim+2];
        for (iDim = 0; iDim < nDim; iDim++)
          Vel_inlet[iDim] = V_inlet[iDim+1];
        
        Vel_normal_inlet_ = 0.0; Vel_inlet_ = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Vel_normal_inlet[iDim] = -Vel_inlet[iDim]*UnitNormal[iDim];
          Vel_normal_inlet_ += Vel_normal_inlet[iDim]*Vel_normal_inlet[iDim];
          Vel_inlet_+= Vel_inlet[iDim]*Vel_inlet[iDim];
        }
        Vel_inlet_ = sqrt(Vel_inlet_);
        Vel_normal_inlet_ = sqrt(Vel_normal_inlet_);
        
        Vel_tangent_inlet_ = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Vel_tangent_inlet[iDim] = Vel_inlet[iDim] - Vel_normal_inlet[iDim];
          Vel_tangent_inlet_ += Vel_tangent_inlet[iDim]*Vel_tangent_inlet[iDim];
        }
        Vel_tangent_inlet_ = sqrt(Vel_tangent_inlet_);
        
        /*--- Mass flow conservation (normal direction) and
         no jump in the tangential velocity ---*/
        
        Vel_normal_outlet_ = (1.0-SecondaryFlow/100.0)*(Rho_inlet*Vel_normal_inlet_)/Rho_outlet;
        
        Vel_outlet_ = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Vel_normal_outlet[iDim] = -Vel_normal_outlet_*UnitNormal[iDim];
          Vel_tangent_outlet[iDim] = Vel_tangent_inlet[iDim];
          Vel_outlet[iDim] = Vel_normal_outlet[iDim] + Vel_tangent_outlet[iDim];
          Vel_outlet_ += Vel_outlet[iDim]*Vel_outlet[iDim];
        }
        Vel_outlet_ = sqrt(Vel_outlet_);
        
        Mach_Outlet = min(Vel_outlet_/SoS_outlet, 1.0);
        
        /*--- Reevaluate the Total Pressure and Total Temperature using the
         Fan Face Mach number and the static values from the jum condition ---*/
        
        Factor = 1.0 + 0.5*Mach_Outlet*Mach_Outlet*Gamma_Minus_One;
        P_Total = P_static * pow(Factor, Gamma/Gamma_Minus_One);
        T_Total = T_static * Factor;
        
        /*--- Flow direction using the velocity direction at the outlet  ---*/
        
        if (Vel_outlet_ != 0.0) {
          for (iDim = 0; iDim < nDim; iDim++) Flow_Dir[iDim] = Vel_outlet[iDim]/Vel_outlet_;
        }
        else {
          for (iDim = 0; iDim < nDim; iDim++) Flow_Dir[iDim] = 0.0;
        }
        
        /*--- Store primitives and set some variables for clarity. ---*/
        
        Density = V_domain[nDim+2];
        Velocity2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = V_domain[iDim+1];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }
        Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
        Pressure    = V_domain[nDim+1];
        H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
        SoundSpeed2 = Gamma*Pressure/Density;
        
        /*--- Compute the acoustic Riemann invariant that is extrapolated
         from the domain interior. ---*/
        
        Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
        for (iDim = 0; iDim < nDim; iDim++)
          Riemann += Velocity[iDim]*UnitNormal[iDim];
        
        /*--- Total speed of sound ---*/
        
        SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
        
        /*--- Dot product of normal and flow direction. This should
         be negative due to outward facing boundary normal convention. ---*/
        
        alpha = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          alpha += UnitNormal[iDim]*Flow_Dir[iDim];
        
        /*--- Coefficients in the quadratic equation for the velocity ---*/
        
        aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
        bb = -1.0*Gamma_Minus_One*alpha*Riemann;
        cc =  0.5*Gamma_Minus_One*Riemann*Riemann - 2.0*SoundSpeed_Total2/Gamma_Minus_One;
        
        /*--- Solve quadratic equation for velocity magnitude. Value must
         be positive, so the choice of root is clear. ---*/
        
        dd = bb*bb - 4.0*aa*cc;
        dd = sqrt(max(0.0, dd));
        Vel_Mag   = (-bb + dd)/(2.0*aa);
        Vel_Mag   = max(0.0, Vel_Mag);
        Velocity2 = Vel_Mag*Vel_Mag;
        
        /*--- Compute speed of sound from total speed of sound eqn. ---*/
        
        SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
        
        /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
        
        Mach2 = min(1.0, Velocity2/SoundSpeed2);
        Velocity2   = Mach2*SoundSpeed2;
        Vel_Mag     = sqrt(Velocity2);
        SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
        
        /*--- Compute new velocity vector at the exit ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
        
        /*--- Static temperature from the speed of sound relation ---*/
        
        Temperature = SoundSpeed2/(Gamma*Gas_Constant);
        
        /*--- Static pressure using isentropic relation at a point ---*/
        
        Pressure = P_Total*pow((Temperature/T_Total), Gamma/Gamma_Minus_One);
        
        /*--- Density at the inlet from the gas law ---*/
        
        Density = Pressure/(Gas_Constant*Temperature);
        
        /*--- Using pressure, density, & velocity, compute the energy ---*/
        
        Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();
        
        /*--- Primitive variables, using the derived quantities ---*/
        
        V_outlet[0] = Temperature;
        for (iDim = 0; iDim < nDim; iDim++)
          V_outlet[iDim+1] = Velocity[iDim];
        V_outlet[nDim+1] = Pressure;
        V_outlet[nDim+2] = Density;
        V_outlet[nDim+3] = Energy + Pressure/Density;
        V_outlet[nDim+4] = sqrt(SoundSpeed2);
        conv_numerics->SetPrimitive(V_domain, V_outlet);
        
      }
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      if (viscous) {
//
//        /*--- Set laminar and eddy viscosity at the infinity ---*/
//
//        if (inlet_surface) {
//          V_inlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
//          V_inlet[nDim+6] = node[iPoint]->GetEddyViscosity();
//        }
//        else {
//          V_outlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
//          V_outlet[nDim+6] = node[iPoint]->GetEddyViscosity();
//        }
//
//        /*--- Set the normal vector and the coordinates ---*/
//
//        visc_numerics->SetNormal(Normal);
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint_Normal]->GetCoord());
//
//        /*--- Primitive variables, and gradient ---*/
//
//        if (inlet_surface) visc_numerics->SetPrimitive(V_domain, V_inlet);
//        else visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
//
//        /*--- Turbulent kinetic energy ---*/
//
//        if (config->GetKind_Turb_Model() == SST)
//          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
//
//        /*--- Compute and update residual ---*/
//
//        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.SubtractBlock(iPoint, Residual);
//
//        /*--- Jacobian contribution for implicit integration ---*/
//
//        if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//
//      }
      
    }
    
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] Flow_Dir;
  
}

void CEulerSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container,
                                CConfig *config, unsigned short val_marker) { }

void CEulerSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

void CEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                        unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;
  
  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
  su2double Volume_nM1, Volume_nP1, TimeStep;
  su2double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Store the physical time step ---*/
  
  TimeStep = config->GetDelta_UnstTimeND();
  
  /*--- Compute the dual time-stepping source term for static meshes ---*/
  
  if (!grid_movement) {
    
    /*--- Loop over all nodes (excluding halos) ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/
      
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                            +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
    
  }
  
  else {
    
    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      /*--- Get indices for nodes i & j plus the face normal ---*/
      
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      Normal = geometry->edge[iEdge]->GetNormal();
      
      /*--- Grid velocities stored at nodes i & j ---*/
      
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      
      /*--- Compute the GCL term by averaging the grid velocities at the
       edge mid-point and dotting with the face normal. ---*/
      
      Residual_GCL = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual_GCL += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      
      /*--- Compute the GCL component of the source term for node i ---*/
      
      U_time_n = node[iPoint]->GetSolution_time_n();
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Compute the GCL component of the source term for node j ---*/
      
      U_time_n = node[jPoint]->GetSolution_time_n();
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      LinSysRes.SubtractBlock(jPoint, Residual);
      
    }
    
    /*---   Loop over the boundary edges ---*/
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        /*--- Get the index for node i plus the boundary face normal ---*/
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        /*--- Grid velocities stored at boundary node i ---*/
        
        GridVel_i = geometry->node[iPoint]->GetGridVel();
        
        /*--- Compute the GCL term by dotting the grid velocity with the face
         normal. The normal is negated to match the boundary convention. ---*/
        
        Residual_GCL = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];
        
        /*--- Compute the GCL component of the source term for node i ---*/
        
        U_time_n = node[iPoint]->GetSolution_time_n();
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;
        LinSysRes.AddBlock(iPoint, Residual);
        
      }
    }
    
    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/
      
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
          + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1/TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (3.0*Volume_nP1)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }
  
}

void CEulerSolver::SetFlow_Displacement(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement,
                                        CConfig *flow_config, CConfig *fea_config, CGeometry **fea_geometry, CSolver ***fea_solution) {
    unsigned short iDim;
    unsigned long iVertex;
    su2double *Coord, VarCoord[3] = {0,0,0};

  #ifndef HAVE_MPI
    unsigned long iPoint_Donor, iPoint;
    unsigned short iMarker;
    su2double *CoordDonor, *DisplacementDonor;

    for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {

      if (flow_config->GetMarker_All_ZoneInterface(iMarker) != 0) {

        for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {

          iPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();

          iPoint_Donor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetDonorPoint();

          Coord = flow_geometry[MESH_0]->node[iPoint]->GetCoord();

          CoordDonor = fea_geometry[MESH_0]->node[iPoint_Donor]->GetCoord();

          /*--- The displacements come from the predicted solution ---*/
          DisplacementDonor = fea_solution[MESH_0][FEA_SOL]->node[iPoint_Donor]->GetSolution_Pred();

          for (iDim = 0; iDim < nDim; iDim++)

            VarCoord[iDim] = (CoordDonor[iDim]+DisplacementDonor[iDim])-Coord[iDim];

          flow_geometry[MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
        }
      }
    }

  #else

    unsigned long nLocalVertexStruct = 0, nLocalVertexFlow = 0;

    unsigned short nMarkerFSI, nMarkerStruct, nMarkerFlow;      // Number of markers on FSI problem, FEA and Flow side
    unsigned short iMarkerFSI, iMarkerStruct, iMarkerFlow;      // Variables for iteration over markers
    unsigned long MaxLocalVertexStruct = 0, MaxLocalVertexFlow = 0;

    unsigned long nBuffer_StructCoord = 0, nBuffer_FlowNewCoord = 0;
    unsigned long nBuffer_DonorIndices = 0, nBuffer_SetIndex = 0;

    unsigned long Point_Flow, Point_Struct;
    long Point_Flow_Rcv, Processor_Flow_Rcv;
    unsigned long Processor_Flow;

    int Marker_Flow = -1, Marker_Struct = -1;

    int iProcessor, nProcessor = 0;


    su2double *Coord_Struct, *Displacement_Struct;

    /*--- Number of markers on the FSI interface ---*/

    nMarkerFSI     = (flow_config->GetMarker_n_ZoneInterface())/2;
    nMarkerStruct  = fea_geometry[MESH_0]->GetnMarker();
    nMarkerFlow    = flow_geometry[MESH_0]->GetnMarker();

    nProcessor = size;

    /*--- Outer loop over the markers on the FSI interface: compute one by one ---*/
    /*--- The tags are always an integer greater than 1: loop from 1 to nMarkerFSI ---*/

    for (iMarkerFSI = 1; iMarkerFSI <= nMarkerFSI; iMarkerFSI++) {

        Marker_Struct = -1;
        Marker_Flow = -1;

        /*--- Initialize pointer buffers inside the loop, so we can delete for each marker. ---*/
        unsigned long Buffer_Send_nVertexStruct[1], *Buffer_Recv_nVertexStruct = NULL;
        unsigned long Buffer_Send_nVertexFlow[1], *Buffer_Recv_nVertexFlow = NULL;

        /*--- The markers on the fluid and structural side are tagged with the same index.
         *--- This is independent of the MPI domain decomposition.
         *--- We need to loop over all markers on structural side and get the number of nodes
         *--- that belong to each FSI marker for each processor ---*/

        /*--- On the structural side ---*/

        for (iMarkerStruct = 0; iMarkerStruct < nMarkerStruct; iMarkerStruct++) {
            /*--- If the tag GetMarker_All_ZoneInterface(iMarkerFEA) equals the index we are looping at ---*/
            if ( fea_config->GetMarker_All_ZoneInterface(iMarkerStruct) == iMarkerFSI ) {
                /*--- We have identified the local index of the FEA marker ---*/
                /*--- Store the number of local points that belong to markFEA on each processor ---*/
                /*--- This includes the halo nodes ---*/
                nLocalVertexStruct = fea_geometry[MESH_0]->GetnVertex(iMarkerStruct);
                /*--- Store the identifier for the structural marker ---*/
                Marker_Struct = iMarkerStruct;
                /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
                break;
            }
            else {
                /*--- If the tag hasn't matched any tag within the FEA markers ---*/
                nLocalVertexStruct = 0;
                Marker_Struct = -1;
            }
        }

        /*--- On the fluid side ---*/

        for (iMarkerFlow = 0; iMarkerFlow < nMarkerFlow; iMarkerFlow++) {
            /*--- If the tag GetMarker_All_ZoneInterface(iMarkerFlow) equals the index we are looping at ---*/
            if ( flow_config->GetMarker_All_ZoneInterface(iMarkerFlow) == iMarkerFSI ) {
                /*--- We have identified the local index of the Flow marker ---*/
                /*--- Store the number of local points that belong to markFlow on each processor ---*/
                /*--- This includes the halo nodes ---*/
                nLocalVertexFlow = flow_geometry[MESH_0]->GetnVertex(iMarkerFlow);
                /*--- Store the identifier for the fluid marker ---*/
                Marker_Flow = iMarkerFlow;
                /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
                break;
            }
            else {
                /*--- If the tag hasn't matched any tag within the Flow markers ---*/
                nLocalVertexFlow = 0;
                Marker_Flow = -1;
            }
        }

        Buffer_Send_nVertexStruct[0] = nLocalVertexStruct;                               // Retrieve total number of vertices on FEA marker
        Buffer_Send_nVertexFlow[0] = nLocalVertexFlow;                               // Retrieve total number of vertices on Flow marker
        if (rank == MASTER_NODE) Buffer_Recv_nVertexStruct = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side
        if (rank == MASTER_NODE) Buffer_Recv_nVertexFlow = new unsigned long[size];  // Allocate memory to receive how many vertices are on each rank on the fluid side

        /*--- We receive MaxLocalVertexFEA as the maximum number of vertices in one single processor on the structural side---*/
        SU2_MPI::Allreduce(&nLocalVertexStruct, &MaxLocalVertexStruct, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        /*--- We receive MaxLocalVertexFlow as the maximum number of vertices in one single processor on the fluid side ---*/
        SU2_MPI::Allreduce(&nLocalVertexFlow, &MaxLocalVertexFlow, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

        /*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the structural side ---*/
        SU2_MPI::Gather(&Buffer_Send_nVertexStruct, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexStruct, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
        /*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the fluid side ---*/
        SU2_MPI::Gather(&Buffer_Send_nVertexFlow, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexFlow, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

        /*--- We will be gathering the structural coordinates into the master node ---*/
        /*--- Then we will distribute them using a scatter operation into the appropriate fluid processor ---*/
        nBuffer_StructCoord = MaxLocalVertexStruct * nDim;
        nBuffer_FlowNewCoord = MaxLocalVertexFlow * nDim;

        /*--- We will be gathering donor index and donor processor (for structure -> donor = flow) ---*/
        /*--- Then we will pass on to the fluid side the index (flow point) to the appropriate processor ---*/
        nBuffer_DonorIndices = 2 * MaxLocalVertexStruct;
        nBuffer_SetIndex = MaxLocalVertexFlow;

        /*--- Send and Recv buffers ---*/

        /*--- Buffers to send and receive the structural coordinates ---*/
        su2double *Buffer_Send_StructCoord = new su2double[nBuffer_StructCoord];
        su2double *Buffer_Recv_StructCoord = NULL;

        /*--- Buffers to send and receive the donor index and processor ---*/
        long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
        long *Buffer_Recv_DonorIndices = NULL;

        /*--- Buffers to send and receive the new fluid coordinates ---*/
        su2double *Buffer_Send_FlowNewCoord = NULL;
        su2double *Buffer_Recv_FlowNewCoord = new su2double[nBuffer_FlowNewCoord];

        /*--- Buffers to send and receive the fluid index ---*/
        long *Buffer_Send_SetIndex = NULL;
        long *Buffer_Recv_SetIndex = new long[nBuffer_SetIndex];

        /*--- Prepare the receive buffers (1st step) and send buffers (2nd step) on the master node only. ---*/

        if (rank == MASTER_NODE) {
            Buffer_Recv_StructCoord  = new su2double[size*nBuffer_StructCoord];
            Buffer_Recv_DonorIndices = new long[size*nBuffer_DonorIndices];
            Buffer_Send_FlowNewCoord = new su2double[size*nBuffer_FlowNewCoord];
            Buffer_Send_SetIndex     = new long[size*nBuffer_SetIndex];
        }

        /*--- On the structural side ---*/

        /*--- If this processor owns the marker we are looping at on the structural side ---*/

        /*--- First we initialize all of the indices and processors to -1 ---*/
        /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
        for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
            Buffer_Send_DonorIndices[iVertex] = -1;

        if (Marker_Struct >= 0) {

            /*--- We have identified the local index of the FEA marker ---*/
            /*--- We loop over all the vertices in that marker and in that particular processor ---*/

            for (iVertex = 0; iVertex < nLocalVertexStruct; iVertex++) {

                Point_Struct = fea_geometry[MESH_0]->vertex[Marker_Struct][iVertex]->GetNode();

                Point_Flow = fea_geometry[MESH_0]->vertex[Marker_Struct][iVertex]->GetDonorPoint();

                Processor_Flow = fea_geometry[MESH_0]->vertex[Marker_Struct][iVertex]->GetDonorProcessor();

                Coord_Struct = fea_geometry[MESH_0]->node[Point_Struct]->GetCoord();

                /*--- The displacements come from the predicted solution ---*/
                Displacement_Struct = fea_solution[MESH_0][FEA_SOL]->node[Point_Struct]->GetSolution_Pred();

                for (iDim = 0; iDim < nDim; iDim++) {
                    Buffer_Send_StructCoord[iVertex*nDim+iDim] = Coord_Struct[iDim] + Displacement_Struct[iDim];
                }
                /*--- If this processor owns the node ---*/
                if (fea_geometry[MESH_0]->node[Point_Struct]->GetDomain()) {
                    Buffer_Send_DonorIndices[2*iVertex]     = Point_Flow;
                    Buffer_Send_DonorIndices[2*iVertex + 1] = Processor_Flow;
                }
                else {
                    /*--- We set the values to be -1 to be able to identify them later as halo nodes ---*/
                    Buffer_Send_DonorIndices[2*iVertex]     = -1;
                    Buffer_Send_DonorIndices[2*iVertex + 1] = -1;
                }

            }
        }

        /*--- Once all the messages have been sent, we gather them all into the MASTER_NODE ---*/
        SU2_MPI::Gather(Buffer_Send_StructCoord, nBuffer_StructCoord, MPI_DOUBLE, Buffer_Recv_StructCoord, nBuffer_StructCoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        SU2_MPI::Gather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);

        /*--- Counter to determine where in the array we have to set the information ---*/
        long *Counter_Processor_Flow = NULL;
        long iProcessor_Struct = 0, iIndex_Struct = 0;
        long iProcessor_Flow = 0, iPoint_Flow = 0, iIndex_Flow = 0;

        /*--- Now we pack the information to send it over to the different processors ---*/

        if (rank == MASTER_NODE) {

            /*--- We set the counter to 0 ---*/
            Counter_Processor_Flow = new long[nProcessor];
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                Counter_Processor_Flow[iProcessor] = 0;
            }

            /*--- First we initialize the index vector to -1 ---*/
            /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
            for (iVertex = 0; iVertex < nProcessor*nBuffer_SetIndex; iVertex++)
                Buffer_Send_SetIndex[iVertex] = -2;

            /*--- As of now we do the loop over the structural points ---*/
            /*--- The number of points for flow and structure does not necessarily have to match ---*/
            /*--- In fact, it's possible that a processor asks for nFlow nodes and there are only ---*/
            /*--- nStruc < nFlow available; this is due to halo nodes ---*/

            /*--- For every processor from which we have received information ---*/
            /*--- (This is, for every processor on the structural side) ---*/
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {

                /*--- This is the initial index on the coordinates buffer for that particular processor on the structural side ---*/
                iProcessor_Struct = iProcessor*nBuffer_StructCoord;
                /*--- This is the initial index on the donor index/processor buffer for that particular processor on the structural side ---*/
                iIndex_Struct = iProcessor*nBuffer_DonorIndices;

                /*--- For every vertex in the information retreived from iProcessor ---*/
                for (iVertex = 0; iVertex < Buffer_Recv_nVertexStruct[iProcessor]; iVertex++) {

                    /*--- The processor and index for the flow are: ---*/
                    Processor_Flow_Rcv = Buffer_Recv_DonorIndices[iIndex_Struct+iVertex*2+1];
                    Point_Flow_Rcv     = Buffer_Recv_DonorIndices[iIndex_Struct+iVertex*2];

                    /*--- Load the buffer at the appropriate position ---*/
                    /*--- This is determined on the fluid side by:
                     *--- Processor_Flow*nBuffer_FlowNewCoord -> Initial position of the processor array (fluid side)
                     *--- +
                     *--- Counter_Processor_Flow*nDim -> Initial position of the nDim array for the particular point on the fluid side
                     *--- +
                     *--- iDim -> Position within the nDim array that corresponds to a point
                     *---
                     *--- While on the structural side is:
                     *--- iProcessor*nBuffer_StructCoord -> Initial position on the processor array (structural side)
                     *--- +
                     *--- iVertex*nDim -> Initial position of the nDim array for the particular point on the structural side
                     */

                    /*--- We check that we are not setting the value for a halo node ---*/
                    if (Point_Flow_Rcv != -1) {
                        iProcessor_Flow = Processor_Flow_Rcv*nBuffer_FlowNewCoord;
                        iIndex_Flow = Processor_Flow_Rcv*nBuffer_SetIndex;
                        iPoint_Flow = Counter_Processor_Flow[Processor_Flow_Rcv]*nDim;

                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Send_FlowNewCoord[iProcessor_Flow + iPoint_Flow + iDim] = Buffer_Recv_StructCoord[iProcessor_Struct + iVertex*nDim + iDim];

                        /*--- We set the fluid index at an appropriate position matching the coordinates ---*/
                        Buffer_Send_SetIndex[iIndex_Flow + Counter_Processor_Flow[Processor_Flow_Rcv]] = Point_Flow_Rcv;

                        Counter_Processor_Flow[Processor_Flow_Rcv]++;
                    }

                }

            }

        }

        /*--- Once all the messages have been prepared, we scatter them all from the MASTER_NODE ---*/
        SU2_MPI::Scatter(Buffer_Send_FlowNewCoord, nBuffer_FlowNewCoord, MPI_DOUBLE, Buffer_Recv_FlowNewCoord, nBuffer_FlowNewCoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        SU2_MPI::Scatter(Buffer_Send_SetIndex, nBuffer_SetIndex, MPI_LONG, Buffer_Recv_SetIndex, nBuffer_SetIndex, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);

        long indexPoint_iVertex, Point_Flow_Check;

        /*--- For the flow marker we are studying ---*/
        if (Marker_Flow >= 0) {

            /*--- We have identified the local index of the Flow marker ---*/
            /*--- We loop over all the vertices in that marker and in that particular processor ---*/

            for (iVertex = 0; iVertex < nLocalVertexFlow; iVertex++) {

                Point_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNode();

                if (flow_geometry[MESH_0]->node[Point_Flow]->GetDomain()) {
                    /*--- Find the index of the point Point_Flow in the buffer Buffer_Recv_SetIndex ---*/
                    indexPoint_iVertex = std::distance(Buffer_Recv_SetIndex, std::find(Buffer_Recv_SetIndex, Buffer_Recv_SetIndex + MaxLocalVertexFlow, Point_Flow));

                    Point_Flow_Check = Buffer_Recv_SetIndex[indexPoint_iVertex];

                    if (Point_Flow_Check < 0) {
                      SU2_MPI::Error("A nonphysical point is being considered for mesh deformation.", CURRENT_FUNCTION);
                    }

                    Coord = flow_geometry[MESH_0]->node[Point_Flow]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                        VarCoord[iDim] = (Buffer_Recv_FlowNewCoord[indexPoint_iVertex*nDim+iDim])-Coord[iDim];

                    flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->SetVarCoord(VarCoord);

                }

            }

        }

        delete [] Buffer_Send_StructCoord;
        delete [] Buffer_Send_DonorIndices;
        delete [] Buffer_Recv_FlowNewCoord;
        delete [] Buffer_Recv_SetIndex;

        if (rank == MASTER_NODE) {
            delete [] Buffer_Recv_nVertexStruct;
            delete [] Buffer_Recv_nVertexFlow;
            delete [] Buffer_Recv_StructCoord;
            delete [] Buffer_Recv_DonorIndices;
            delete [] Buffer_Send_FlowNewCoord;
            delete [] Buffer_Send_SetIndex;
            delete [] Counter_Processor_Flow;
        }

    }

  #endif

    flow_grid_movement->SetVolume_Deformation(flow_geometry[MESH_0], flow_config, true);

}

void CEulerSolver::SetFlow_Displacement_Int(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement,
                                        CConfig *flow_config, CConfig *fea_config, CGeometry **fea_geometry, CSolver ***fea_solution) {
    unsigned short iMarker, iDim, iDonor, nDonor;
    unsigned long iVertex;
    su2double VarCoord[3];

    unsigned long iPoint_Donor;
    su2double *DisplacementDonor, *DisplacementDonor_Prev, coeff;

    for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {

      if (flow_config->GetMarker_All_ZoneInterface(iMarker) != 0) {

        for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {

          for (iDim = 0; iDim < nDim; iDim++)
            VarCoord[iDim]=0.0;

          nDonor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetnDonorPoints();

          for (iDonor = 0; iDonor < nDonor; iDonor++) {
            iPoint_Donor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetInterpDonorPoint(iDonor);
            coeff = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetDonorCoeff(iDonor);

            /*--- The displacements come from the predicted solution ---*/
            DisplacementDonor = fea_solution[MESH_0][FEA_SOL]->node[iPoint_Donor]->GetSolution_Pred();

            DisplacementDonor_Prev = fea_solution[MESH_0][FEA_SOL]->node[iPoint_Donor]->GetSolution_Pred_Old();

            for (iDim = 0; iDim < nDim; iDim++)

              VarCoord[iDim] += (DisplacementDonor[iDim] - DisplacementDonor_Prev[iDim])*coeff;
          }

          flow_geometry[MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
        }
      }
    }
    flow_grid_movement->SetVolume_Deformation(flow_geometry[MESH_0], flow_config, true);

}

void CEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  unsigned short turb_model = config->GetKind_Turb_Model();
  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine;
  bool grid_movement  = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool static_fsi = ((config->GetUnsteady_Simulation() == STEADY) &&
                     (config->GetFSI_Simulation()));
  bool steady_restart = config->GetSteadyRestart();
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  string UnstExt, text_line;
  ifstream restart_file;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  string restart_filename = config->GetSolution_FlowFileName();

  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Skip coordinates ---*/
  
  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/

  if (dual_time || time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      node[iPoint_Local]->SetSolution(Solution);
      iPoint_Global_Local++;

      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/

      if (grid_movement && val_update_geo) {

        /*--- First, remove any variables for the turbulence model that
         appear in the restart file before the grid velocities. ---*/

        if (turb_model == SA || turb_model == SA_NEG) {
          index++;
        } else if (turb_model == SST) {
          index+=2;
        }

        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
        /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/

        su2double GridVel[3] = {0.0,0.0,0.0};
        if (!steady_restart) {

          /*--- Rewind the index to retrieve the Coords. ---*/
          index = counter*Restart_Vars[1];
          for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim]; }

          /*--- Move the index forward to get the grid velocities. ---*/
          index = counter*Restart_Vars[1] + skipVars + nVar;
          for (iDim = 0; iDim < nDim; iDim++) { GridVel[iDim] = Restart_Data[index+iDim]; }
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
          geometry[MESH_0]->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
        }
      }
        
      if (static_fsi && val_update_geo) {
       /*--- Rewind the index to retrieve the Coords. ---*/
        index = counter*Restart_Vars[1];
        for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim];}

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
        }
      }

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We alo call the preprocessing routine
   on the fine level in order to have all necessary quantities updated,
   especially if this is a turbulent simulation (eddy viscosity). ---*/

  solver[MESH_0][FLOW_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][FLOW_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(Solution);
    }
    solver[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  }

  /*--- Update the geometry for flows on dynamic meshes ---*/

  if (grid_movement && val_update_geo) {

    /*--- Communicate the new coordinates and grid velocities at the halos ---*/

    geometry[MESH_0]->Set_MPI_Coord(config);
    geometry[MESH_0]->Set_MPI_GridVel(config);

    /*--- Recompute the edges and dual mesh control volumes in the
     domain and on the boundaries. ---*/

    geometry[MESH_0]->SetCoord_CG();
    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);

    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/

    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config, geometry[iMeshFine], UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config, geometry[iMeshFine],UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
      geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMeshFine], config);
      }
    }

  /*--- Update the geometry for flows on static FSI problems with moving meshes ---*/

  if (static_fsi && val_update_geo) {

    /*--- Communicate the new coordinates and grid velocities at the halos ---*/

    geometry[MESH_0]->Set_MPI_Coord(config);

    /*--- Recompute the edges and  dual mesh control volumes in the
     domain and on the boundaries. ---*/

    geometry[MESH_0]->SetCoord_CG();
    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);

    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/

    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config, geometry[iMeshFine], UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config, geometry[iMeshFine],UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
    }
  }
  
  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

void CEulerSolver::SetFreeStream_Solution(CConfig *config) {

  unsigned long iPoint;
  unsigned short iDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->SetSolution(0, Density_Inf);
    for (iDim = 0; iDim < nDim; iDim++) {
      node[iPoint]->SetSolution(iDim+1, Density_Inf*Velocity_Inf[iDim]);
    }
    node[iPoint]->SetSolution(nVar-1, Density_Inf*Energy_Inf);
  }
}

void CEulerSolver::SetFreeStream_TurboSolution(CConfig *config) {

  unsigned long iPoint;
  unsigned short iDim;
  unsigned short iZone  =  config->GetiZone();
  su2double *turboVelocity, *cartVelocity, *turboNormal;

  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Mach             = config->GetMach();
  su2double SoundSpeed;

  turboVelocity   = new su2double[nDim];
  cartVelocity    = new su2double[nDim];

  turboNormal     = config->GetFreeStreamTurboNormal();

  FluidModel->SetTDState_Prho(Pressure_Inf, Density_Inf);
  SoundSpeed = FluidModel->GetSoundSpeed();

  /*--- Compute the Free Stream velocity, using the Mach number ---*/
  turboVelocity[0] = cos(Alpha)*Mach*SoundSpeed;
  turboVelocity[1] = sin(Alpha)*Mach*SoundSpeed;


  if (nDim == 3) {
    turboVelocity[2] = 0.0;
  }

  ComputeBackVelocity(turboVelocity, turboNormal, cartVelocity, INFLOW, config->GetKind_TurboMachinery(iZone));

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->SetSolution(0, Density_Inf);
    for (iDim = 0; iDim < nDim; iDim++) {
      node[iPoint]->SetSolution(iDim+1, Density_Inf*cartVelocity[iDim]);
    }
    node[iPoint]->SetSolution(nVar-1, Density_Inf*Energy_Inf);

    node[iPoint]->SetPrimVar(FluidModel);
    node[iPoint]->SetSecondaryVar(FluidModel);
  }

  delete [] turboVelocity;
  delete [] cartVelocity;


}



void CEulerSolver::PreprocessAverage(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, iMarkerTP, iSpan;
  su2double Pressure = 0.0, Density = 0.0, *Velocity = NULL, *TurboVelocity,
      Area, TotalArea, TotalAreaPressure, TotalAreaDensity, *TotalAreaVelocity, *UnitNormal, *TurboNormal;
  string Marker_Tag, Monitoring_Tag;
  unsigned short  iZone     = config->GetiZone();
  su2double  *AverageTurboNormal, VelSq;

  /*-- Variables declaration and allocation ---*/
  Velocity           = new su2double[nDim];
  UnitNormal         = new su2double[nDim];
  TurboNormal        = new su2double[nDim];
  TurboVelocity      = new su2double[nDim];
  TotalAreaVelocity  = new su2double[nDim];

#ifdef HAVE_MPI
  su2double MyTotalAreaDensity, MyTotalAreaPressure;
  su2double *MyTotalAreaVelocity = NULL;
#endif

  for (iSpan= 0; iSpan < nSpanWiseSections; iSpan++){


    for (iDim=0; iDim<nDim; iDim++) {
      TotalAreaVelocity[iDim] = 0.0;
    }

    TotalAreaPressure = 0.0;
    TotalAreaDensity  = 0.0;

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){

            /*--- Retrieve Old Solution ---*/

            /*--- Loop over the vertices to sum all the quantities pithc-wise ---*/
            for (iVertex = 0; iVertex < geometry->GetnVertexSpan(iMarker,iSpan); iVertex++) {
              iPoint = geometry->turbovertex[iMarker][iSpan][iVertex]->GetNode();
              if (geometry->node[iPoint]->GetDomain()){
                /*--- Compute the integral fluxes for the boundaries ---*/

                Pressure = node[iPoint]->GetPressure();
                Density = node[iPoint]->GetDensity();

                /*--- Normal vector for this vertex (negate for outward convention) ---*/
                geometry->turbovertex[iMarker][iSpan][iVertex]->GetNormal(UnitNormal);
                geometry->turbovertex[iMarker][iSpan][iVertex]->GetTurboNormal(TurboNormal);
                Area = geometry->turbovertex[iMarker][iSpan][iVertex]->GetArea();

                VelSq = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                  Velocity[iDim] = node[iPoint]->GetPrimitive(iDim+1);
                  VelSq += Velocity[iDim]*Velocity[iDim];
                }

                ComputeTurboVelocity(Velocity, TurboNormal , TurboVelocity, marker_flag, config->GetKind_TurboMachinery(iZone));

                /*--- Compute different integral quantities for the boundary of interest ---*/

                TotalAreaPressure += Area*Pressure;
                TotalAreaDensity  += Area*Density;
                for (iDim = 0; iDim < nDim; iDim++)
                  TotalAreaVelocity[iDim] += Area*Velocity[iDim];
              }
            }
          }
        }
      }
    }


#ifdef HAVE_MPI

    /*--- Add information using all the nodes ---*/

    MyTotalAreaDensity   = TotalAreaDensity;          TotalAreaDensity     = 0;
    MyTotalAreaPressure  = TotalAreaPressure;         TotalAreaPressure    = 0;

    SU2_MPI::Allreduce(&MyTotalAreaDensity, &TotalAreaDensity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalAreaPressure, &TotalAreaPressure, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    MyTotalAreaVelocity    = new su2double[nDim];


    for (iDim = 0; iDim < nDim; iDim++) {
      MyTotalAreaVelocity[iDim]    = TotalAreaVelocity[iDim];
      TotalAreaVelocity[iDim]      = 0.0;
    }

    SU2_MPI::Allreduce(MyTotalAreaVelocity, TotalAreaVelocity, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    delete [] MyTotalAreaVelocity;


#endif

    /*--- initialize spanwise average quantities ---*/


    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){

            TotalArea           = geometry->GetSpanArea(iMarker,iSpan);
            AverageTurboNormal 	= geometry->GetAverageTurboNormal(iMarker,iSpan);

            /*--- Compute the averaged value for the boundary of interest for the span of interest ---*/

            AverageDensity[iMarker][iSpan]           = TotalAreaDensity / TotalArea;
            AveragePressure[iMarker][iSpan]          = TotalAreaPressure / TotalArea;
            for (iDim = 0; iDim < nDim; iDim++)
              AverageVelocity[iMarker][iSpan][iDim]  = TotalAreaVelocity[iDim] / TotalArea;

            /* --- compute static averaged quantities ---*/
            ComputeTurboVelocity(AverageVelocity[iMarker][iSpan], AverageTurboNormal , AverageTurboVelocity[iMarker][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));

            OldAverageDensity[iMarker][iSpan]               = AverageDensity[iMarker][iSpan];
            OldAveragePressure[iMarker][iSpan]              = AveragePressure[iMarker][iSpan];
            for(iDim = 0; iDim < nDim;iDim++)
              OldAverageTurboVelocity[iMarker][iSpan][iDim] = AverageTurboVelocity[iMarker][iSpan][iDim];

          }
        }
      }
    }
  }

  /*--- initialize 1D average quantities ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
      if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){

          AverageTurboNormal 	= geometry->GetAverageTurboNormal(iMarker,nSpanWiseSections);

          /*--- Compute the averaged value for the boundary of interest for the span of interest ---*/

          AverageDensity[iMarker][nSpanWiseSections]          = AverageDensity[iMarker][nSpanWiseSections/2];
          AveragePressure[iMarker][nSpanWiseSections]         = AveragePressure[iMarker][nSpanWiseSections/2];
          for (iDim = 0; iDim < nDim; iDim++)
            AverageVelocity[iMarker][nSpanWiseSections][iDim] = AverageVelocity[iMarker][nSpanWiseSections/2][iDim];

          /* --- compute static averaged quantities ---*/
          ComputeTurboVelocity(AverageVelocity[iMarker][nSpanWiseSections], AverageTurboNormal , AverageTurboVelocity[iMarker][nSpanWiseSections], marker_flag, config->GetKind_TurboMachinery(iZone));

          OldAverageDensity[iMarker][nSpanWiseSections]               = AverageDensity[iMarker][nSpanWiseSections];
          OldAveragePressure[iMarker][nSpanWiseSections]              = AveragePressure[iMarker][nSpanWiseSections];
          for(iDim = 0; iDim < nDim;iDim++)
            OldAverageTurboVelocity[iMarker][nSpanWiseSections][iDim] = AverageTurboVelocity[iMarker][nSpanWiseSections][iDim];

        }
      }
    }
  }


  /*--- Free locally allocated memory ---*/
  delete [] Velocity;
  delete [] UnitNormal;
  delete [] TurboNormal;
  delete [] TurboVelocity;
  delete [] TotalAreaVelocity;

}


void CEulerSolver::TurboAverageProcess(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag) {

  unsigned long iVertex, iPoint, nVert;
  unsigned short iDim, iVar, iMarker, iMarkerTP, iSpan, jSpan;
  unsigned short average_process = config->GetKind_AverageProcess();
  unsigned short performance_average_process = config->GetKind_PerformanceAverageProcess();
  su2double Pressure = 0.0, Density = 0.0, Enthalpy = 0.0,  *Velocity = NULL, *TurboVelocity,
      Area, TotalArea, Radius1, Radius2, Vt2, TotalAreaPressure, TotalAreaDensity, *TotalAreaVelocity, *UnitNormal, *TurboNormal,
      TotalMassPressure, TotalMassDensity, *TotalMassVelocity;
  string Marker_Tag, Monitoring_Tag;
  su2double val_init_pressure;
  unsigned short  iZone     = config->GetiZone();
  su2double TotalDensity, TotalPressure, *TotalVelocity, *AverageTurboNormal, *TotalFluxes;
  su2double TotalNu, TotalOmega, TotalKine, TotalMassNu, TotalMassOmega, TotalMassKine, TotalAreaNu, TotalAreaOmega, TotalAreaKine;
  su2double Nu, Kine, Omega;
  su2double MachTest, soundSpeed;
  bool turbulent = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool spalart_allmaras = (config->GetKind_Turb_Model() == SA);
  bool menter_sst       = (config->GetKind_Turb_Model() == SST);

  /*-- Variables declaration and allocation ---*/
  Velocity            = new su2double[nDim];
  UnitNormal          = new su2double[nDim];
  TurboNormal         = new su2double[nDim];
  TurboVelocity       = new su2double[nDim];
  TotalVelocity       = new su2double[nDim];
  TotalAreaVelocity   = new su2double[nDim];
  TotalMassVelocity   = new su2double[nDim];
  TotalFluxes         = new su2double[nVar];

  su2double avgDensity, *avgVelocity, avgPressure, avgKine, avgOmega, avgNu, avgAreaDensity, *avgAreaVelocity, avgAreaPressure,
  avgAreaKine, avgAreaOmega, avgAreaNu, avgMassDensity, *avgMassVelocity, avgMassPressure, avgMassKine, avgMassOmega, avgMassNu,
  avgMixDensity, *avgMixVelocity, *avgMixTurboVelocity, avgMixPressure, avgMixKine, avgMixOmega, avgMixNu;

  avgVelocity         = new su2double[nDim];
  avgAreaVelocity     = new su2double[nDim];
  avgMassVelocity     = new su2double[nDim];
  avgMixVelocity      = new su2double[nDim];
  avgMixTurboVelocity = new su2double[nDim];


#ifdef HAVE_MPI
  su2double MyTotalDensity, MyTotalPressure, MyTotalAreaDensity, MyTotalAreaPressure, MyTotalMassPressure, MyTotalMassDensity, *MyTotalFluxes = NULL;
  su2double *MyTotalVelocity = NULL, *MyTotalAreaVelocity = NULL, *MyTotalMassVelocity = NULL;
  su2double MyTotalNu, MyTotalKine, MyTotalOmega, MyTotalAreaNu, MyTotalAreaKine, MyTotalAreaOmega, MyTotalMassNu, MyTotalMassKine, MyTotalMassOmega;
#endif


  for (iSpan= 0; iSpan < nSpanWiseSections + 1; iSpan++){

    /*--- Forces initialization for contenitors ---*/
    for (iVar=0;iVar<nVar;iVar++)
      TotalFluxes[iVar]= 0.0;
    for (iDim=0; iDim<nDim; iDim++) {
      TotalVelocity[iDim]     = 0.0;
      TotalAreaVelocity[iDim] = 0.0;
      TotalMassVelocity[iDim] = 0.0;
    }

    TotalDensity      = 0.0;
    TotalPressure     = 0.0;
    TotalAreaPressure = 0.0;
    TotalAreaDensity  = 0.0;
    TotalMassPressure = 0.0;
    TotalMassDensity  = 0.0;
    TotalNu           = 0.0;
    TotalOmega        = 0.0;
    TotalKine         = 0.0;
    TotalMassNu       = 0.0;
    TotalMassOmega    = 0.0;
    TotalMassKine     = 0.0;
    TotalAreaNu       = 0.0;
    TotalAreaOmega    = 0.0;
    TotalAreaKine     = 0.0;

    Nu    = 0.0;
    Omega = 0.0;
    Kine  = 0.0;


    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){

            /*--- Retrieve Old Solution ---*/

            /*--- Loop over the vertices to sum all the quantities pithc-wise ---*/
            if(iSpan < nSpanWiseSections){
              for (iVertex = 0; iVertex < geometry->GetnVertexSpan(iMarker,iSpan); iVertex++) {
                iPoint = geometry->turbovertex[iMarker][iSpan][iVertex]->GetNode();

                /*--- Compute the integral fluxes for the boundaries ---*/
                Pressure = node[iPoint]->GetPressure();
                Density  = node[iPoint]->GetDensity();
                Enthalpy = node[iPoint]->GetEnthalpy();

                /*--- Normal vector for this vertex (negate for outward convention) ---*/
                geometry->turbovertex[iMarker][iSpan][iVertex]->GetNormal(UnitNormal);
                geometry->turbovertex[iMarker][iSpan][iVertex]->GetTurboNormal(TurboNormal);
                Area = geometry->turbovertex[iMarker][iSpan][iVertex]->GetArea();
                su2double VelNormal = 0.0, VelSq = 0.0;

                for (iDim = 0; iDim < nDim; iDim++) {
                  Velocity[iDim] = node[iPoint]->GetPrimitive(iDim+1);
                  VelNormal += UnitNormal[iDim]*Velocity[iDim];
                  VelSq += Velocity[iDim]*Velocity[iDim];
                }

                ComputeTurboVelocity(Velocity, TurboNormal , TurboVelocity, marker_flag, config->GetKind_TurboMachinery(iZone));

                /*--- Compute different integral quantities for the boundary of interest ---*/

                TotalDensity          += Density;
                TotalPressure         += Pressure;
                for (iDim = 0; iDim < nDim; iDim++)
                  TotalVelocity[iDim] += Velocity[iDim];

                TotalAreaPressure         += Area*Pressure;
                TotalAreaDensity          += Area*Density;
                for (iDim = 0; iDim < nDim; iDim++)
                  TotalAreaVelocity[iDim] += Area*Velocity[iDim];

                TotalMassPressure         += Area*(Density*TurboVelocity[0] )*Pressure;
                TotalMassDensity          += Area*(Density*TurboVelocity[0] )*Density;
                for (iDim = 0; iDim < nDim; iDim++)
                  TotalMassVelocity[iDim] += Area*(Density*TurboVelocity[0] )*Velocity[iDim];

                TotalFluxes[0]      += Area*(Density*TurboVelocity[0]);
                TotalFluxes[1]      += Area*(Density*TurboVelocity[0]*TurboVelocity[0] + Pressure);
                for (iDim = 2; iDim < nDim+1; iDim++)
                  TotalFluxes[iDim] += Area*(Density*TurboVelocity[0]*TurboVelocity[iDim -1]);
                TotalFluxes[nDim+1] += Area*(Density*TurboVelocity[0]*Enthalpy);


                /*--- Compute turbulent integral quantities for the boundary of interest ---*/

                if(turbulent){
                  if(menter_sst){
                    Kine = solver[TURB_SOL]->node[iPoint]->GetSolution(0);
                    Omega = solver[TURB_SOL]->node[iPoint]->GetSolution(1);
                  }
                  if(spalart_allmaras){
                    Nu = solver[TURB_SOL]->node[iPoint]->GetSolution(0);
                  }

                  TotalKine   += Kine;
                  TotalOmega  += Omega;
                  TotalNu     += Nu;

                  TotalAreaKine    += Area*Kine;
                  TotalAreaOmega   += Area*Omega;
                  TotalAreaNu      += Area*Nu;

                  TotalMassKine    += Area*(Density*TurboVelocity[0] )*Kine;
                  TotalMassOmega   += Area*(Density*TurboVelocity[0] )*Omega;
                  TotalMassNu      += Area*(Density*TurboVelocity[0] )*Nu;


                }
              }
            }
            else{
              for (jSpan= 0; jSpan < nSpanWiseSections; jSpan++){
                for (iVertex = 0; iVertex < geometry->GetnVertexSpan(iMarker,jSpan); iVertex++) {
                  iPoint = geometry->turbovertex[iMarker][jSpan][iVertex]->GetNode();

                  /*--- Compute the integral fluxes for the boundaries ---*/
                  Pressure = node[iPoint]->GetPressure();
                  Density  = node[iPoint]->GetDensity();
                  Enthalpy = node[iPoint]->GetEnthalpy();

                  /*--- Normal vector for this vertex (negate for outward convention) ---*/
                  geometry->turbovertex[iMarker][jSpan][iVertex]->GetNormal(UnitNormal);
                  geometry->turbovertex[iMarker][jSpan][iVertex]->GetTurboNormal(TurboNormal);
                  Area = geometry->turbovertex[iMarker][jSpan][iVertex]->GetArea();
                  su2double VelNormal = 0.0, VelSq = 0.0;

                  for (iDim = 0; iDim < nDim; iDim++) {
                    Velocity[iDim] = node[iPoint]->GetPrimitive(iDim+1);
                    VelNormal += UnitNormal[iDim]*Velocity[iDim];
                    VelSq += Velocity[iDim]*Velocity[iDim];
                  }

                  ComputeTurboVelocity(Velocity, TurboNormal , TurboVelocity, marker_flag, config->GetKind_TurboMachinery(iZone));

                  /*--- Compute different integral quantities for the boundary of interest ---*/

                  TotalDensity          += Density;
                  TotalPressure         += Pressure;
                  for (iDim = 0; iDim < nDim; iDim++)
                    TotalVelocity[iDim] += Velocity[iDim];

                  TotalAreaPressure         += Area*Pressure;
                  TotalAreaDensity          += Area*Density;
                  for (iDim = 0; iDim < nDim; iDim++)
                    TotalAreaVelocity[iDim] += Area*Velocity[iDim];

                  TotalMassPressure         += Area*(Density*TurboVelocity[0] )*Pressure;
                  TotalMassDensity          += Area*(Density*TurboVelocity[0] )*Density;
                  for (iDim = 0; iDim < nDim; iDim++)
                    TotalMassVelocity[iDim] += Area*(Density*TurboVelocity[0] )*Velocity[iDim];

                  TotalFluxes[0]      += Area*(Density*TurboVelocity[0]);
                  TotalFluxes[1]      += Area*(Density*TurboVelocity[0]*TurboVelocity[0] + Pressure);
                  for (iDim = 2; iDim < nDim+1; iDim++)
                    TotalFluxes[iDim] += Area*(Density*TurboVelocity[0]*TurboVelocity[iDim -1]);
                  TotalFluxes[nDim+1] += Area*(Density*TurboVelocity[0]*Enthalpy);


                  /*--- Compute turbulent integral quantities for the boundary of interest ---*/

                  if(turbulent){
                    if(menter_sst){
                      Kine  = solver[TURB_SOL]->node[iPoint]->GetSolution(0);
                      Omega = solver[TURB_SOL]->node[iPoint]->GetSolution(1);
                    }
                    if(spalart_allmaras){
                      Nu    = solver[TURB_SOL]->node[iPoint]->GetSolution(0);
                    }

                    TotalKine   += Kine;
                    TotalOmega  += Omega;
                    TotalNu     += Nu;

                    TotalAreaKine   += Area*Kine;
                    TotalAreaOmega  += Area*Omega;
                    TotalAreaNu     += Area*Nu;

                    TotalMassKine    += Area*(Density*TurboVelocity[0] )*Kine;
                    TotalMassOmega   += Area*(Density*TurboVelocity[0] )*Omega;
                    TotalMassNu      += Area*(Density*TurboVelocity[0] )*Nu;

                  }
                }
              }
            }
          }
        }
      }
    }

#ifdef HAVE_MPI

    /*--- Add information using all the nodes ---*/

    MyTotalDensity       = TotalDensity;              TotalDensity         = 0;
    MyTotalPressure      = TotalPressure;             TotalPressure        = 0;
    MyTotalAreaDensity   = TotalAreaDensity;          TotalAreaDensity     = 0;
    MyTotalAreaPressure  = TotalAreaPressure;         TotalAreaPressure    = 0;
    MyTotalMassDensity   = TotalMassDensity;          TotalMassDensity     = 0;
    MyTotalMassPressure  = TotalMassPressure;         TotalMassPressure    = 0;
    MyTotalNu            = TotalNu;                   TotalNu              = 0;
    MyTotalKine          = TotalKine;                 TotalKine            = 0;
    MyTotalOmega         = TotalOmega;                TotalOmega           = 0;
    MyTotalAreaNu        = TotalAreaNu;               TotalAreaNu          = 0;
    MyTotalAreaKine      = TotalAreaKine;             TotalAreaKine        = 0;
    MyTotalAreaOmega     = TotalAreaOmega;            TotalAreaOmega       = 0;
    MyTotalMassNu        = TotalMassNu;               TotalMassNu          = 0;
    MyTotalMassKine      = TotalMassKine;             TotalMassKine        = 0;
    MyTotalMassOmega     = TotalMassOmega;            TotalMassOmega       = 0;




    SU2_MPI::Allreduce(&MyTotalDensity, &TotalDensity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalPressure, &TotalPressure, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalAreaDensity, &TotalAreaDensity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalAreaPressure, &TotalAreaPressure, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalMassDensity, &TotalMassDensity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalMassPressure, &TotalMassPressure, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalNu, &TotalNu, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalKine, &TotalKine, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalOmega, &TotalOmega, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalAreaNu, &TotalAreaNu, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalAreaKine, &TotalAreaKine, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalAreaOmega, &TotalAreaOmega, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalMassNu, &TotalMassNu, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalMassKine, &TotalMassKine, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalMassOmega, &TotalMassOmega, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    MyTotalFluxes          = new su2double[nVar];
    MyTotalVelocity        = new su2double[nDim];
    MyTotalAreaVelocity    = new su2double[nDim];
    MyTotalMassVelocity    = new su2double[nDim];

    for (iVar = 0; iVar < nVar; iVar++) {
      MyTotalFluxes[iVar]  = TotalFluxes[iVar];
      TotalFluxes[iVar]    = 0.0;
    }

    for (iDim = 0; iDim < nDim; iDim++) {
      MyTotalVelocity[iDim]        = TotalVelocity[iDim];
      TotalVelocity[iDim]          = 0.0;
      MyTotalAreaVelocity[iDim]    = TotalAreaVelocity[iDim];
      TotalAreaVelocity[iDim]      = 0.0;
      MyTotalMassVelocity[iDim]    = TotalMassVelocity[iDim];
      TotalMassVelocity[iDim]      = 0.0;
    }

    SU2_MPI::Allreduce(MyTotalFluxes, TotalFluxes, nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MyTotalVelocity, TotalVelocity, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MyTotalAreaVelocity, TotalAreaVelocity, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MyTotalMassVelocity, TotalMassVelocity, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    delete [] MyTotalFluxes; delete [] MyTotalVelocity; delete [] MyTotalAreaVelocity; delete [] MyTotalMassVelocity;


#endif


    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){

            TotalArea           = geometry->GetSpanArea(iMarker,iSpan);
            AverageTurboNormal 	= geometry->GetAverageTurboNormal(iMarker,iSpan);
            nVert               = geometry->GetnTotVertexSpan(iMarker,iSpan);



            /*--- compute normal Mach number as a check for massflow average and mixedout average ---*/
            FluidModel->SetTDState_Prho(TotalAreaPressure/TotalArea, TotalAreaDensity / TotalArea);
            soundSpeed = FluidModel->GetSoundSpeed();
            MachTest   = TotalFluxes[0]/(TotalAreaDensity*soundSpeed);

            /*--- Compute the averaged value for the boundary of interest for the span of interest ---*/

            /*--- compute algebraic average ---*/
            avgDensity        = TotalDensity / nVert;
            avgPressure       = TotalPressure / nVert;
            for (iDim = 0; iDim < nDim; iDim++) avgVelocity[iDim] = TotalVelocity[iDim] / nVert;
            avgKine           = TotalKine/nVert;
            avgOmega          = TotalOmega/nVert;
            avgNu             = TotalNu/nVert;

            /*--- compute area average ---*/
            avgAreaDensity     = TotalAreaDensity / TotalArea;
            avgAreaPressure    = TotalAreaPressure / TotalArea;
            for (iDim = 0; iDim < nDim; iDim++) avgAreaVelocity[iDim] = TotalAreaVelocity[iDim] / TotalArea;
            avgAreaKine        = TotalAreaKine / TotalArea;
            avgAreaOmega       = TotalAreaOmega / TotalArea;
            avgAreaNu          = TotalAreaNu / TotalArea;

            /*--- compute mass-flow average ---*/
            if (abs(MachTest)< config->GetAverageMachLimit()) {
              avgMassDensity   = avgAreaDensity;
              avgMassPressure  = avgAreaPressure;
              for (iDim = 0; iDim < nDim; iDim++) avgMassVelocity[iDim] = avgAreaVelocity[iDim];
              avgMassKine      = avgAreaKine;
              avgMassOmega     = avgAreaOmega;
              avgMassNu        = avgAreaNu;
            }else{
              avgMassDensity     = TotalMassDensity / TotalFluxes[0];
              avgMassPressure    = TotalMassPressure / TotalFluxes[0];
              for (iDim = 0; iDim < nDim; iDim++) avgMassVelocity[iDim] = TotalMassVelocity[iDim] / TotalFluxes[0];
              avgMassKine        = TotalMassKine / TotalFluxes[0];
              avgMassOmega       = TotalMassOmega / TotalFluxes[0];
              avgMassNu          = TotalMassNu / TotalFluxes[0];
            }
            /*--- compute mixed-out average ---*/
            for (iVar = 0; iVar<nVar; iVar++){
              AverageFlux[iMarker][iSpan][iVar]   = TotalFluxes[iVar]/TotalArea;
              SpanTotalFlux[iMarker][iSpan][iVar] = TotalFluxes[iVar];
            }
            val_init_pressure = OldAveragePressure[iMarker][iSpan];

            if (abs(MachTest)< config->GetAverageMachLimit()) {
              avgMixDensity    = avgAreaDensity;
              avgMixPressure   = avgAreaPressure;
              for (iDim = 0; iDim < nDim; iDim++)
                avgMixVelocity[iDim] = avgAreaVelocity[iDim];
              ComputeTurboVelocity(avgMixVelocity, AverageTurboNormal , avgMixTurboVelocity, marker_flag, config->GetKind_TurboMachinery(iZone));
              avgMixKine       = avgAreaKine;
              avgMixOmega      = avgAreaOmega;
              avgMixNu         = avgAreaNu;
            }else {
              MixedOut_Average (config, val_init_pressure, AverageFlux[iMarker][iSpan], AverageTurboNormal, avgMixPressure, avgMixDensity);
              avgMixTurboVelocity[0]         = ( AverageFlux[iMarker][iSpan][1] - avgMixPressure) / AverageFlux[iMarker][iSpan][0];
              for (iDim = 2; iDim < nDim +1;iDim++)
                avgMixTurboVelocity[iDim-1]  = AverageFlux[iMarker][iSpan][iDim] / AverageFlux[iMarker][iSpan][0];

              if (avgMixDensity!= avgMixDensity || avgMixPressure!= avgMixPressure || avgMixPressure < 0.0 || avgMixDensity < 0.0 ){
                val_init_pressure = avgAreaPressure;
                MixedOut_Average (config, val_init_pressure, AverageFlux[iMarker][iSpan], AverageTurboNormal, avgMixPressure, avgMixDensity);
                avgMixTurboVelocity[0]          = ( AverageFlux[iMarker][iSpan][1] - avgMixPressure) / AverageFlux[iMarker][iSpan][0];
                for (iDim = 2; iDim < nDim +1;iDim++)
                  avgMixTurboVelocity[iDim-1]   = AverageFlux[iMarker][iSpan][iDim] / AverageFlux[iMarker][iSpan][0];
              }
              avgMixKine       = avgMassKine;
              avgMixOmega      = avgMassOmega;
              avgMixNu         = avgMassNu;
            }

            /*--- Store averaged value for the selected average method ---*/
            switch(average_process){
            case ALGEBRAIC:
              AverageDensity[iMarker][iSpan]  = avgDensity;
              AveragePressure[iMarker][iSpan] = avgPressure;
              ComputeTurboVelocity(avgVelocity, AverageTurboNormal , AverageTurboVelocity[iMarker][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
              AverageKine[iMarker][iSpan]     = avgKine;
              AverageOmega[iMarker][iSpan]    = avgOmega;
              AverageNu[iMarker][iSpan]       = avgNu;
              break;

            case AREA:
              AverageDensity[iMarker][iSpan]  = avgAreaDensity;
              AveragePressure[iMarker][iSpan] = avgAreaPressure;
              ComputeTurboVelocity(avgAreaVelocity, AverageTurboNormal , AverageTurboVelocity[iMarker][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
              AverageKine[iMarker][iSpan]     = avgAreaKine;
              AverageOmega[iMarker][iSpan]    = avgAreaOmega;
              AverageNu[iMarker][iSpan]       = avgAreaNu;
              break;

            case MASSFLUX:
              AverageDensity[iMarker][iSpan]  = avgMassDensity;
              AveragePressure[iMarker][iSpan] = avgMassPressure;
              ComputeTurboVelocity(avgAreaVelocity, AverageTurboNormal , AverageTurboVelocity[iMarker][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
              AverageKine[iMarker][iSpan]     = avgMassKine;
              AverageOmega[iMarker][iSpan]    = avgMassOmega;
              AverageNu[iMarker][iSpan]       = avgMassNu;
              break;

            case MIXEDOUT:
              AverageDensity[iMarker][iSpan] = avgMixDensity;
              AveragePressure[iMarker][iSpan] = avgMixPressure;
              for (iDim = 0; iDim < nDim; iDim++) AverageTurboVelocity[iMarker][iSpan][iDim] = avgMixTurboVelocity[iDim];
              AverageKine[iMarker][iSpan]     = avgMixKine;
              AverageOmega[iMarker][iSpan]    = avgMixOmega;
              AverageNu[iMarker][iSpan]       = avgMixNu;
              break;

            default:
              SU2_MPI::Error(" Invalid AVERAGE PROCESS input!", CURRENT_FUNCTION);
              break;
            }

            /* --- check if averaged quantities are correct otherwise reset the old quantities ---*/
            if (AverageDensity[iMarker][iSpan]!= AverageDensity[iMarker][iSpan] || AveragePressure[iMarker][iSpan]!= AveragePressure[iMarker][iSpan]){
              cout<<"nan in mixing process routine for iSpan: " << iSpan<< " in marker " << config->GetMarker_All_TagBound(iMarker)<< endl;
              AverageDensity[iMarker][iSpan]               = OldAverageDensity[iMarker][iSpan];
              AveragePressure[iMarker][iSpan]              = OldAveragePressure[iMarker][iSpan];
              for(iDim = 0; iDim < nDim;iDim++)
                AverageTurboVelocity[iMarker][iSpan][iDim] = OldAverageTurboVelocity[iMarker][iSpan][iDim];
            }


            if (AverageDensity[iMarker][iSpan] < 0.0 || AveragePressure[iMarker][iSpan] < 0.0){
              cout << " negative density or pressure in mixing process routine for iSpan: " << iSpan<< " in marker " << config->GetMarker_All_TagBound(iMarker)<< endl;
              AverageDensity[iMarker][iSpan]               = OldAverageDensity[iMarker][iSpan];
              AveragePressure[iMarker][iSpan]              = OldAveragePressure[iMarker][iSpan];
              for(iDim = 0; iDim < nDim;iDim++)
                AverageTurboVelocity[iMarker][iSpan][iDim] = OldAverageTurboVelocity[iMarker][iSpan][iDim];
            }

            /* --- update old average solution ---*/
            OldAverageDensity[iMarker][iSpan]               = AverageDensity[iMarker][iSpan];
            OldAveragePressure[iMarker][iSpan]              = AveragePressure[iMarker][iSpan];
            for(iDim = 0; iDim < nDim;iDim++)
              OldAverageTurboVelocity[iMarker][iSpan][iDim] = AverageTurboVelocity[iMarker][iSpan][iDim];

            /*--- to avoid back flow ---*/
            if (AverageTurboVelocity[iMarker][iSpan][0] < 0.0){
              AverageTurboVelocity[iMarker][iSpan][0]       = soundSpeed*config->GetAverageMachLimit();
            }

            /*--- compute cartesian average Velocity ---*/
            ComputeBackVelocity(AverageTurboVelocity[iMarker][iSpan], AverageTurboNormal , AverageVelocity[iMarker][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));


            /*--- Store averaged performance value for the selected average method ---*/
            switch(performance_average_process){
            case ALGEBRAIC:
              if(marker_flag == INFLOW){
                DensityIn[iMarkerTP - 1][iSpan]   = avgDensity;
                PressureIn[iMarkerTP - 1][iSpan]  = avgPressure;
                ComputeTurboVelocity(avgVelocity, AverageTurboNormal , TurboVelocityIn[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                KineIn[iMarkerTP - 1][iSpan]      = avgKine;
                OmegaIn[iMarkerTP - 1][iSpan]     = avgOmega;
                NuIn[iMarkerTP - 1][iSpan]        = avgNu;
              }
              else{
                DensityOut[iMarkerTP - 1][iSpan]  = avgDensity;
                PressureOut[iMarkerTP - 1][iSpan] = avgPressure;
                ComputeTurboVelocity(avgVelocity, AverageTurboNormal , TurboVelocityOut[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                KineOut[iMarkerTP - 1][iSpan]     = avgKine;
                OmegaOut[iMarkerTP - 1][iSpan]    = avgOmega;
                NuOut[iMarkerTP - 1][iSpan]       = avgNu;
              }

              break;
            case AREA:
              if(marker_flag == INFLOW){
                DensityIn[iMarkerTP - 1][iSpan]   = avgAreaDensity;
                PressureIn[iMarkerTP - 1][iSpan]  = avgAreaPressure;
                ComputeTurboVelocity(avgAreaVelocity, AverageTurboNormal , TurboVelocityIn[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                KineIn[iMarkerTP - 1][iSpan]      = avgAreaKine;
                OmegaIn[iMarkerTP - 1][iSpan]     = avgAreaOmega;
                NuIn[iMarkerTP - 1][iSpan]        = avgAreaNu;
              }
              else{
                DensityOut[iMarkerTP - 1][iSpan]  = avgAreaDensity;
                PressureOut[iMarkerTP - 1][iSpan] = avgAreaPressure;
                ComputeTurboVelocity(avgAreaVelocity, AverageTurboNormal , TurboVelocityOut[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                KineOut[iMarkerTP - 1][iSpan]     = avgAreaKine;
                OmegaOut[iMarkerTP - 1][iSpan]    = avgAreaOmega;
                NuOut[iMarkerTP - 1][iSpan]       = avgAreaNu/TotalArea;
              }
              break;

            case MASSFLUX:
              if(marker_flag == INFLOW){
                DensityIn[iMarkerTP - 1][iSpan]   = avgMassDensity;
                PressureIn[iMarkerTP - 1][iSpan]  = avgMassPressure;
                ComputeTurboVelocity(avgMassVelocity, AverageTurboNormal , TurboVelocityIn[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                KineIn[iMarkerTP - 1][iSpan]      = avgMassKine;
                OmegaIn[iMarkerTP - 1][iSpan]     = avgMassOmega;
                NuIn[iMarkerTP - 1][iSpan]        = avgMassNu;
              }
              else{
                DensityOut[iMarkerTP - 1][iSpan]  = avgMassDensity;
                PressureOut[iMarkerTP - 1][iSpan] = avgMassPressure;
                ComputeTurboVelocity(avgMassVelocity, AverageTurboNormal , TurboVelocityOut[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                KineOut[iMarkerTP - 1][iSpan]     = avgMassKine;
                OmegaOut[iMarkerTP - 1][iSpan]    = avgMassOmega;
                NuOut[iMarkerTP - 1][iSpan]       = avgMassNu;
              }

              break;

            case MIXEDOUT:
              if (marker_flag == INFLOW){
                DensityIn[iMarkerTP - 1][iSpan]  = avgMixDensity;
                PressureIn[iMarkerTP - 1][iSpan] = avgMixPressure;
                for (iDim = 0; iDim < nDim; iDim++) TurboVelocityIn[iMarkerTP -1][iSpan][iDim] = avgMixTurboVelocity[iDim];
                KineIn[iMarkerTP - 1][iSpan]     = avgMixKine;
                OmegaIn[iMarkerTP - 1][iSpan]    = avgMixOmega;
                NuIn[iMarkerTP - 1][iSpan]       = avgMixNu;
              }
              else{
                DensityOut[iMarkerTP - 1][iSpan]  = avgMixDensity;
                PressureOut[iMarkerTP - 1][iSpan] = avgMixPressure;
                for (iDim = 0; iDim < nDim; iDim++) TurboVelocityOut[iMarkerTP -1][iSpan][iDim] = avgMixTurboVelocity[iDim];
                KineOut[iMarkerTP - 1][iSpan]     = avgMixKine;
                OmegaOut[iMarkerTP - 1][iSpan]    = avgMixOmega;
                NuOut[iMarkerTP - 1][iSpan]       = avgMixNu;
              }
              break;

            default:
              SU2_MPI::Error(" Invalid MIXING_PROCESS input!", CURRENT_FUNCTION);
              break;
            }
          }
        }
      }
    }
  }



  /*--- Compute Outlet Static Pressure if Radial equilibrium is imposed ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
      if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
          Marker_Tag         = config->GetMarker_All_TagBound(iMarker);
          if(config->GetBoolGiles() || config->GetBoolRiemann()){
            if(config->GetBoolRiemann()){
              if(config->GetKind_Data_Riemann(Marker_Tag) == RADIAL_EQUILIBRIUM){
                RadialEquilibriumPressure[iMarker][nSpanWiseSections/2] = config->GetRiemann_Var1(Marker_Tag)/config->GetPressure_Ref();
                for (iSpan= nSpanWiseSections/2; iSpan < nSpanWiseSections-1; iSpan++){
                  Radius2    = geometry->GetTurboRadius(iMarker,iSpan+1);
                  Radius1    = geometry->GetTurboRadius(iMarker,iSpan);
                  Vt2        = AverageTurboVelocity[iMarker][iSpan +1][1]*AverageTurboVelocity[iMarker][iSpan +1][1];
                  RadialEquilibriumPressure[iMarker][iSpan +1] =  RadialEquilibriumPressure[iMarker][iSpan] + AverageDensity[iMarker][iSpan +1]*Vt2/Radius2*(Radius2 - Radius1);
                }
                for (iSpan= nSpanWiseSections/2; iSpan > 0; iSpan--){
                  Radius2    = geometry->GetTurboRadius(iMarker,iSpan);
                  Radius1    = geometry->GetTurboRadius(iMarker,iSpan-1);
                  Vt2        = AverageTurboVelocity[iMarker][iSpan - 1][1]*AverageTurboVelocity[iMarker][iSpan - 1][1];
                  RadialEquilibriumPressure[iMarker][iSpan -1] =  RadialEquilibriumPressure[iMarker][iSpan] - AverageDensity[iMarker][iSpan -1]*Vt2/Radius2*(Radius2 - Radius1);
                }
              }
            }
            else{
              if(config->GetKind_Data_Giles(Marker_Tag) == RADIAL_EQUILIBRIUM){
                RadialEquilibriumPressure[iMarker][nSpanWiseSections/2] = config->GetGiles_Var1(Marker_Tag)/config->GetPressure_Ref();
                for (iSpan= nSpanWiseSections/2; iSpan < nSpanWiseSections-1; iSpan++){
                  Radius2    = geometry->GetTurboRadius(iMarker,iSpan+1);
                  Radius1    = geometry->GetTurboRadius(iMarker,iSpan);
                  Vt2        = AverageTurboVelocity[iMarker][iSpan +1][1]*AverageTurboVelocity[iMarker][iSpan +1][1];
                  RadialEquilibriumPressure[iMarker][iSpan +1] =  RadialEquilibriumPressure[iMarker][iSpan] + AverageDensity[iMarker][iSpan +1]*Vt2/Radius2*(Radius2 - Radius1);
                }
                for (iSpan= nSpanWiseSections/2; iSpan > 0; iSpan--){
                  Radius2    = geometry->GetTurboRadius(iMarker,iSpan);
                  Radius1    = geometry->GetTurboRadius(iMarker,iSpan-1);
                  Vt2        = AverageTurboVelocity[iMarker][iSpan -1][1]*AverageTurboVelocity[iMarker][iSpan - 1][1];
                  RadialEquilibriumPressure[iMarker][iSpan -1] =  RadialEquilibriumPressure[iMarker][iSpan] - AverageDensity[iMarker][iSpan -1]*Vt2/Radius2*(Radius2 - Radius1);
                }
              }
            }
          }
        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Velocity;
  delete [] UnitNormal;
  delete [] TurboNormal;
  delete [] TurboVelocity;
  delete [] TotalVelocity;
  delete [] TotalAreaVelocity;
  delete [] TotalFluxes;
  delete [] TotalMassVelocity;
  delete [] avgVelocity;
  delete [] avgAreaVelocity;
  delete [] avgMassVelocity;
  delete [] avgMixVelocity;
  delete [] avgMixTurboVelocity;

}

void CEulerSolver::MixedOut_Average (CConfig *config, su2double val_init_pressure, su2double *val_Averaged_Flux, su2double *val_normal,
    su2double& pressure_mix, su2double& density_mix) {


  su2double dx, f, df, resdl = 1.0E+05;
  unsigned short iter = 0, iDim;
  su2double relax_factor = config->GetMixedout_Coeff(0);
  su2double toll = config->GetMixedout_Coeff(1);
  unsigned short maxiter = SU2_TYPE::Int(config->GetMixedout_Coeff(2));
  su2double dhdP, dhdrho, enthalpy_mix, velsq, *vel;

  vel = new su2double[nDim];

  pressure_mix = val_init_pressure;

  /*--- Newton-Raphson's method with central difference formula ---*/

  while ( iter <= maxiter ) {

    density_mix = val_Averaged_Flux[0]*val_Averaged_Flux[0]/(val_Averaged_Flux[1] - pressure_mix);
    FluidModel->SetTDState_Prho(pressure_mix, density_mix);
    enthalpy_mix = FluidModel->GetStaticEnergy() + (pressure_mix)/(density_mix);

    FluidModel->ComputeDerivativeNRBC_Prho(pressure_mix, density_mix);
    dhdP   = FluidModel->GetdhdP_rho();
    dhdrho = FluidModel->Getdhdrho_P();

    vel[0]  = (val_Averaged_Flux[1] - pressure_mix) / val_Averaged_Flux[0];
    for (iDim = 1; iDim < nDim; iDim++) {
      vel[iDim]  = val_Averaged_Flux[iDim+1] / val_Averaged_Flux[0];
    }

    velsq = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      velsq += vel[iDim]*vel[iDim];
    }
    f = val_Averaged_Flux[nDim+1] - val_Averaged_Flux[0]*(enthalpy_mix + velsq/2);
    df = -val_Averaged_Flux[0]*(dhdP - 1/density_mix) - dhdrho*density_mix*density_mix/val_Averaged_Flux[0];
    dx = -f/df;
    resdl = dx/val_init_pressure;
    pressure_mix += relax_factor*dx;

    iter += 1;
    if ( abs(resdl) <= toll ) {
      break;
    }

  }

  density_mix = val_Averaged_Flux[0]*val_Averaged_Flux[0]/(val_Averaged_Flux[1] - pressure_mix);


  delete [] vel;

}

void CEulerSolver::GatherInOutAverageValues(CConfig *config, CGeometry *geometry){

  unsigned short iMarker, iMarkerTP;
  unsigned short iSpan;
  int markerTP;
  su2double     densityIn, pressureIn, normalVelocityIn, tangVelocityIn, radialVelocityIn;
  su2double     densityOut, pressureOut, normalVelocityOut, tangVelocityOut, radialVelocityOut;
  su2double     kineIn, omegaIn, nuIn, kineOut, omegaOut, nuOut;
  //TODO (turbo) implement interpolation so that Inflow and Outflow spanwise section can be different

  for (iSpan= 0; iSpan < nSpanWiseSections + 1 ; iSpan++) {
    
#ifdef HAVE_MPI
    unsigned short i, n1, n2, n1t,n2t;
    su2double *TurbPerfIn= NULL,*TurbPerfOut= NULL;
    su2double *TotTurbPerfIn = NULL,*TotTurbPerfOut = NULL;
    int *TotMarkerTP = NULL;

    n1          = 8;
    n2          = 8;
    n1t         = n1*size;
    n2t         = n2*size;
    TurbPerfIn  = new su2double[n1];
    TurbPerfOut = new su2double[n2];

    for (i=0;i<n1;i++)
      TurbPerfIn[i]    = -1.0;
    for (i=0;i<n2;i++)
      TurbPerfOut[i]   = -1.0;
#endif


    densityIn            = -1.0;
    pressureIn           = -1.0;
    normalVelocityIn     = -1.0;
    tangVelocityIn       = -1.0;
    radialVelocityIn     = -1.0;
    densityOut           = -1.0;
    pressureOut          = -1.0;
    normalVelocityOut    = -1.0;
    tangVelocityOut      = -1.0;
    radialVelocityOut    = -1.0;
    kineIn               = -1.0;
    omegaIn              = -1.0;
    nuIn                 = -1.0;
    kineOut              = -1.0;
    omegaOut             = -1.0;
    nuOut                = -1.0;


    markerTP             = -1;

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      for (iMarkerTP = 1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == INFLOW){
            markerTP            = iMarkerTP;
            densityIn           = DensityIn[iMarkerTP -1][iSpan];
            pressureIn          = PressureIn[iMarkerTP -1][iSpan];
            normalVelocityIn    = TurboVelocityIn[iMarkerTP -1][iSpan][0];
            tangVelocityIn      = TurboVelocityIn[iMarkerTP -1][iSpan][1];
            if (nDim ==3){
              radialVelocityIn  = TurboVelocityIn[iMarkerTP -1][iSpan][2];
            }
            kineIn              = KineIn[iMarkerTP -1][iSpan];
            omegaIn             = OmegaIn[iMarkerTP -1][iSpan];
            nuIn                = NuIn[iMarkerTP -1][iSpan];


#ifdef HAVE_MPI
            TurbPerfIn[0]  = densityIn;
            TurbPerfIn[1]  = pressureIn;
            TurbPerfIn[2]  = normalVelocityIn;
            TurbPerfIn[3]  = tangVelocityIn;
            TurbPerfIn[4]  = radialVelocityIn;
            TurbPerfIn[5]  = kineIn;
            TurbPerfIn[6]  = omegaIn;
            TurbPerfIn[7]  = nuIn;
#endif
          }

          /*--- retrieve outlet information ---*/
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == OUTFLOW){
            densityOut           = DensityOut[iMarkerTP -1][iSpan];
            pressureOut          = PressureOut[iMarkerTP -1][iSpan];
            normalVelocityOut    = TurboVelocityOut[iMarkerTP -1][iSpan][0];
            tangVelocityOut      = TurboVelocityOut[iMarkerTP -1][iSpan][1];
            if (nDim ==3){
              radialVelocityOut  = TurboVelocityOut[iMarkerTP -1][iSpan][2];
            }
            kineOut              = KineOut[iMarkerTP -1][iSpan];
            omegaOut             = OmegaOut[iMarkerTP -1][iSpan];
            nuOut                = NuOut[iMarkerTP -1][iSpan];



#ifdef HAVE_MPI
            TurbPerfOut[0]  = densityOut;
            TurbPerfOut[1]  = pressureOut;
            TurbPerfOut[2]  = normalVelocityOut;
            TurbPerfOut[3]  = tangVelocityOut;
            TurbPerfOut[4]  = radialVelocityOut;
            TurbPerfOut[5]  = kineOut;
            TurbPerfOut[6]  = omegaOut;
            TurbPerfOut[7]  = nuOut;
#endif
          }
        }
      }
    }

#ifdef HAVE_MPI
    if (rank == MASTER_NODE){
      TotTurbPerfIn       = new su2double[n1t];
      TotTurbPerfOut      = new su2double[n2t];
      for (i=0;i<n1t;i++)
        TotTurbPerfIn[i]  = -1.0;
      for (i=0;i<n2t;i++)
        TotTurbPerfOut[i] = -1.0;
      TotMarkerTP = new int[size];
      for(i=0; i<size; i++){
        TotMarkerTP[i]    = -1;
      }
    }
    SU2_MPI::Gather(TurbPerfIn, n1, MPI_DOUBLE, TotTurbPerfIn, n1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(TurbPerfOut, n2, MPI_DOUBLE,TotTurbPerfOut, n2, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(&markerTP, 1, MPI_INT,TotMarkerTP, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    if (rank == MASTER_NODE){
      delete [] TurbPerfIn, delete [] TurbPerfOut;
    }

    if (rank == MASTER_NODE){
      for (i=0;i<size;i++){
        if(TotTurbPerfIn[n1*i] > 0.0){
          densityIn        = 0.0;
          densityIn        = TotTurbPerfIn[n1*i];
          pressureIn       = 0.0;
          pressureIn       = TotTurbPerfIn[n1*i+1];
          normalVelocityIn = 0.0;
          normalVelocityIn = TotTurbPerfIn[n1*i+2];
          tangVelocityIn   = 0.0;
          tangVelocityIn   = TotTurbPerfIn[n1*i+3];
          radialVelocityIn = 0.0;
          radialVelocityIn = TotTurbPerfIn[n1*i+4];
          kineIn           = 0.0;
          kineIn           = TotTurbPerfIn[n1*i+5];
          omegaIn          = 0.0;
          omegaIn          = TotTurbPerfIn[n1*i+6];
          nuIn             = 0.0;
          nuIn             = TotTurbPerfIn[n1*i+7];

          markerTP         = -1;
          markerTP         = TotMarkerTP[i];
        }

        if(TotTurbPerfOut[n2*i] > 0.0){
          densityOut        = 0.0;
          densityOut        = TotTurbPerfOut[n1*i];
          pressureOut       = 0.0;
          pressureOut       = TotTurbPerfOut[n1*i+1];
          normalVelocityOut = 0.0;
          normalVelocityOut = TotTurbPerfOut[n1*i+2];
          tangVelocityOut   = 0.0;
          tangVelocityOut   = TotTurbPerfOut[n1*i+3];
          radialVelocityOut = 0.0;
          radialVelocityOut = TotTurbPerfOut[n1*i+4];
          kineOut           = 0.0;
          kineOut           = TotTurbPerfOut[n1*i+5];
          omegaOut          = 0.0;
          omegaOut          = TotTurbPerfOut[n1*i+6];
          nuOut             = 0.0;
          nuOut             = TotTurbPerfOut[n1*i+7];
        }
      }

      delete [] TotTurbPerfIn, delete [] TotTurbPerfOut; delete [] TotMarkerTP;
    }

#endif

    if (rank == MASTER_NODE && markerTP > -1){
      /*----Quantities needed for computing the turbomachinery performance -----*/
      DensityIn[markerTP -1][iSpan]              = densityIn;
      PressureIn[markerTP -1][iSpan]             = pressureIn;
      TurboVelocityIn[markerTP -1][iSpan][0]     = normalVelocityIn;
      TurboVelocityIn[markerTP -1][iSpan][1]     = tangVelocityIn;
      if (nDim == 3)
        TurboVelocityIn[markerTP -1][iSpan][2]   = radialVelocityIn;
      KineIn[markerTP -1][iSpan]                 = kineIn;
      OmegaIn[markerTP -1][iSpan]                = omegaIn;
      NuIn[markerTP -1][iSpan]                   = nuIn;

      DensityOut[markerTP -1][iSpan]             = densityOut;
      PressureOut[markerTP -1][iSpan]            = pressureOut;
      TurboVelocityOut[markerTP -1][iSpan][0]    = normalVelocityOut;
      TurboVelocityOut[markerTP -1][iSpan][1]    = tangVelocityOut;
      if (nDim == 3)
        TurboVelocityOut[markerTP -1][iSpan][2]  = radialVelocityOut;
      KineOut[markerTP -1][iSpan]                = kineOut;
      OmegaOut[markerTP -1][iSpan]               = omegaOut;
      NuOut[markerTP -1][iSpan]                  = nuOut;
    }
  }
}

CNSSolver::CNSSolver(void) : CEulerSolver() {
  
  /*--- Basic array initialization ---*/
  
  CD_Visc = NULL; CL_Visc = NULL; CSF_Visc = NULL; CEff_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CoPx_Visc = NULL;   CoPy_Visc = NULL;   CoPz_Visc = NULL;

  ForceViscous = NULL; MomentViscous = NULL; CSkinFriction = NULL;
  
  /*--- Surface based array initialization ---*/
  
  Surface_CL_Visc = NULL; Surface_CD_Visc = NULL; Surface_CSF_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  Surface_HF_Visc = NULL; Surface_MaxHF_Visc = NULL;
  
  /*--- Rotorcraft simulation array initialization ---*/
  
  CMerit_Visc = NULL; CT_Visc = NULL; CQ_Visc = NULL;

  /*--- Inlet Variables ---*/
  Inlet_Ttotal = NULL;
  Inlet_Ptotal = NULL;
  Inlet_FlowDir = NULL;
  
  SlidingState      = NULL;
  SlidingStateNodes = NULL;
  
}

CNSSolver::CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CEulerSolver() {
  
  unsigned long iPoint, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  su2double Density, Velocity2, Pressure, Temperature, StaticEnergy;
  ifstream restart_file;
  unsigned short nZone = geometry->GetnZone();
  bool restart    = (config->GetRestart() || config->GetRestart_Flow());
  int Unst_RestartIter;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  bool roe_turkel = (config->GetKind_Upwind_Flow() == TURKEL);
  bool low_mach_prec = config->Low_Mach_Preconditioning();
  
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  string filename_ = config->GetSolution_FlowFileName();

  unsigned short direct_diff = config->GetDirectDiff();
  bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS));

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone);

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/

    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Read and store the restart metadata. ---*/

    Read_SU2_Restart_Metadata(geometry, config, false, filename_);
    
  }

  /*--- Array initialization ---*/
  
  CD_Visc = NULL; CL_Visc = NULL; CSF_Visc = NULL; CEff_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CoPx_Visc = NULL;   CoPy_Visc = NULL;   CoPz_Visc = NULL;

  Surface_CL_Visc = NULL; Surface_CD_Visc = NULL; Surface_CSF_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  Surface_HF_Visc = NULL; Surface_MaxHF_Visc = NULL;
  
  CMerit_Visc = NULL;      CT_Visc = NULL;      CQ_Visc = NULL;
  MaxHF_Visc = NULL; ForceViscous = NULL; MomentViscous = NULL;
  CSkinFriction = NULL;    Cauchy_Serie = NULL; HF_Visc = NULL;
  
  /*--- Initialize quantities for the average process for internal flow ---*/

  AverageVelocity 				    = NULL;
  AverageTurboVelocity 				= NULL;
  ExtAverageTurboVelocity 			= NULL;
  AverageFlux 						= NULL;
  SpanTotalFlux 					= NULL;
  AveragePressure  					= NULL;
  RadialEquilibriumPressure         = NULL;
  ExtAveragePressure  				= NULL;
  AverageDensity   					= NULL;
  ExtAverageDensity   				= NULL;
  AverageNu                         = NULL;
  AverageKine                       = NULL;
  AverageOmega                      = NULL;
  ExtAverageNu                      = NULL;
  ExtAverageKine                    = NULL;
  ExtAverageOmega                   = NULL;


  /*--- Initialize primitive quantities for turboperformace ---*/

  DensityIn                     = NULL;
  PressureIn                    = NULL;
  TurboVelocityIn               = NULL;
  DensityOut                    = NULL;
  PressureOut                   = NULL;
  TurboVelocityOut              = NULL;


  /*--- Initialize quantities for Giles BC ---*/

  CkInflow                      = NULL;
  CkOutflow1			= NULL;
  CkOutflow2			= NULL;



  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure
   Compressible flow, primitive variables (T, vx, vy, vz, P, rho, h, c, lamMu, EddyMu, ThCond, Cp).
   ---*/
  
  nDim = geometry->GetnDim();
  
  nVar = nDim+2;
  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
  nSecondaryVar = 8; nSecondaryVarGrad = 2;

  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nPrimVarGrad;
  
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
 
  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex = new unsigned long[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];
 
  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  
  SetNondimensionalization(geometry, config, iMesh);
  
  /*--- Allocate the node variables ---*/
  node = new CVariable*[nPoint];
  
  /*--- Define some auxiliar vector related with the residual ---*/
  
  Residual      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Residual_i    = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
  
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Define some auxiliary vectors related to the primitive solution ---*/
  
  Primitive   = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the Secondary solution ---*/
  
  Secondary   = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar]   = 0.0;
  Secondary_i = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_i[iVar] = 0.0;
  Secondary_j = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_j[iVar] = 0.0;

  /*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }
  
  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/
  
  if (roe_turkel || low_mach_prec) {
    LowMach_Precontioner = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new su2double[nVar];
  }
  
  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  }
  
  else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Define some auxiliary vectors for computing flow variable
   gradients by least squares, S matrix := inv(R)*traspose(inv(R)),
   c vector := transpose(WA)*(Wb) ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    Cvector = new su2double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  CharacPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }
  
  /*--- Store the value of the primitive variables + 2 turb variables at the boundaries,
   used for IO with a donor cell ---*/
  
  DonorPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      if (rans) {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar+2];
        for (iVar = 0; iVar < nPrimVar + 2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
      else {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
        for (iVar = 0; iVar < nPrimVar ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
    }
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorGlobalIndex[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Delta P at the Actuator Disk ---*/
  
  ActDisk_DeltaP = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ActDisk_DeltaP[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      ActDisk_DeltaP[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Delta T at the Actuator Disk ---*/
  
  ActDisk_DeltaT = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ActDisk_DeltaT[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      ActDisk_DeltaT[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Total Pressure at the inlet BC ---*/
  
  Inlet_Ttotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
      Inlet_Ttotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        Inlet_Ttotal[iMarker][iVertex] = 0;
      }
    } else {
      Inlet_Ttotal[iMarker] = NULL;
    }
  }
  
  /*--- Store the value of the Total Temperature at the inlet BC ---*/
  
  Inlet_Ptotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
      Inlet_Ptotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        Inlet_Ptotal[iMarker][iVertex] = 0;
      }
    } else {
      Inlet_Ptotal[iMarker] = NULL;
    }
  }
  
  /*--- Store the value of the Flow direction at the inlet BC ---*/
  
  Inlet_FlowDir = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
      Inlet_FlowDir[iMarker] = new su2double* [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        Inlet_FlowDir[iMarker][iVertex] = new su2double [nDim];
        for (iDim = 0; iDim < nDim; iDim++) {
          Inlet_FlowDir[iMarker][iVertex][iDim] = 0;
        }
      }
    } else {
      Inlet_FlowDir[iMarker] = NULL;
    }
  }

  /*--- Set up inlet profiles, if necessary ---*/

  SetInlet(config);

  /*--- Inviscid force definition and coefficient in all the markers ---*/
  
  CPressure = new su2double* [nMarker];
  CPressureTarget = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new su2double [geometry->nVertex[iMarker]];
    CPressureTarget[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressure[iMarker][iVertex] = 0.0;
      CPressureTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Heat flux in all the markers ---*/
  
  HeatFlux = new su2double* [nMarker];
  HeatFluxTarget = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    HeatFlux[iMarker] = new su2double [geometry->nVertex[iMarker]];
    HeatFluxTarget[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      HeatFlux[iMarker][iVertex] = 0.0;
      HeatFluxTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Y plus in all the markers ---*/
  
  YPlus = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    YPlus[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      YPlus[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Skin friction in all the markers ---*/
  
  CSkinFriction = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CSkinFriction[iMarker] = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      CSkinFriction[iMarker][iDim] = new su2double[geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        CSkinFriction[iMarker][iDim][iVertex] = 0.0;
      }
    }
  }
  
  /*--- Non dimensional coefficients ---*/
  
  ForceInviscid  = new su2double[3];
  MomentInviscid = new su2double[3];
  CD_Inv         = new su2double[nMarker];
  CL_Inv         = new su2double[nMarker];
  CSF_Inv        = new su2double[nMarker];
  CEff_Inv       = new su2double[nMarker];
  CFx_Inv        = new su2double[nMarker];
  CFy_Inv        = new su2double[nMarker];
  CFz_Inv        = new su2double[nMarker];
  CMx_Inv        = new su2double[nMarker];
  CMy_Inv        = new su2double[nMarker];
  CMz_Inv        = new su2double[nMarker];
  CoPx_Inv        = new su2double[nMarker];
  CoPy_Inv        = new su2double[nMarker];
  CoPz_Inv        = new su2double[nMarker];

  ForceMomentum  = new su2double[3];
  MomentMomentum = new su2double[3];
  CD_Mnt         = new su2double[nMarker];
  CL_Mnt         = new su2double[nMarker];
  CSF_Mnt        = new su2double[nMarker];
  CEff_Mnt       = new su2double[nMarker];
  CFx_Mnt        = new su2double[nMarker];
  CFy_Mnt        = new su2double[nMarker];
  CFz_Mnt        = new su2double[nMarker];
  CMx_Mnt        = new su2double[nMarker];
  CMy_Mnt        = new su2double[nMarker];
  CMz_Mnt        = new su2double[nMarker];
  CoPx_Mnt        = new su2double[nMarker];
  CoPy_Mnt        = new su2double[nMarker];
  CoPz_Mnt        = new su2double[nMarker];

  ForceViscous     = new su2double[3];
  MomentViscous    = new su2double[3];
  CD_Visc          = new su2double[nMarker];
  CL_Visc          = new su2double[nMarker];
  CSF_Visc         = new su2double[nMarker];
  CEff_Visc        = new su2double[nMarker];
  CFx_Visc         = new su2double[nMarker];
  CFy_Visc         = new su2double[nMarker];
  CFz_Visc         = new su2double[nMarker];
  CMx_Visc         = new su2double[nMarker];
  CMy_Visc         = new su2double[nMarker];
  CMz_Visc         = new su2double[nMarker];
  CoPx_Visc         = new su2double[nMarker];
  CoPy_Visc         = new su2double[nMarker];
  CoPz_Visc         = new su2double[nMarker];

  Surface_CL_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Inv     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];
 
  Surface_CL_Mnt         = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Mnt         = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Mnt        = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF         = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff           = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz            = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Visc     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Visc       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_HF_Visc         = new su2double[config->GetnMarker_Monitoring()];
  Surface_MaxHF_Visc      = new su2double[config->GetnMarker_Monitoring()];
  
  /*--- Rotational coefficients ---*/
  
  CMerit_Inv = new su2double[nMarker];
  CT_Inv     = new su2double[nMarker];
  CQ_Inv     = new su2double[nMarker];
  
  CMerit_Mnt = new su2double[nMarker];
  CT_Mnt     = new su2double[nMarker];
  CQ_Mnt     = new su2double[nMarker];

  CMerit_Visc = new su2double[nMarker];
  CT_Visc     = new su2double[nMarker];
  CQ_Visc     = new su2double[nMarker];
  
  /*--- Heat based coefficients ---*/
  
  HF_Visc    = new su2double[nMarker];
  MaxHF_Visc = new su2double[nMarker];
  
  /*--- Supersonic coefficients ---*/
  
  CEquivArea_Inv   = new su2double[nMarker];
  CNearFieldOF_Inv = new su2double[nMarker];
  
  /*--- Engine simulation ---*/
  
  Inflow_MassFlow     = new su2double[nMarker];
  Inflow_Pressure     = new su2double[nMarker];
  Inflow_Mach         = new su2double[nMarker];
  Inflow_Area         = new su2double[nMarker];
  
  Exhaust_MassFlow    = new su2double[nMarker];
  Exhaust_Pressure    = new su2double[nMarker];
  Exhaust_Temperature = new su2double[nMarker];
  Exhaust_Area        = new su2double[nMarker];
  
  /*--- Init total coefficients ---*/
  
  Total_CD         = 0.0;    Total_CL           = 0.0;    Total_CSF          = 0.0;
  Total_CMx        = 0.0;    Total_CMy          = 0.0;    Total_CMz          = 0.0;
  Total_CoPx        = 0.0;    Total_CoPy          = 0.0;    Total_CoPz          = 0.0;
  Total_CEff       = 0.0;    Total_CEquivArea   = 0.0;    Total_CNearFieldOF = 0.0;
  Total_CFx        = 0.0;    Total_CFy          = 0.0;    Total_CFz          = 0.0;
  Total_CT         = 0.0;    Total_CQ           = 0.0;    Total_CMerit       = 0.0;
  Total_MaxHeat    = 0.0;   Total_Heat         = 0.0;    Total_ComboObj     = 0.0;
  Total_CpDiff     = 0.0;   Total_HeatFluxDiff = 0.0;    Total_BCThrust_Prev = 0.0;
  Total_NetCThrust = 0.0;   Total_NetCThrust_Prev = 0.0; Total_CL_Prev = 0.0;
  Total_Power      = 0.0;   AoA_Prev           = 0.0;    Total_CD_Prev      = 0.0;
  Total_CMx_Prev   = 0.0;   Total_CMy_Prev      = 0.0;   Total_CMz_Prev      = 0.0;
  Total_AeroCD     = 0.0;   Total_SolidCD     = 0.0;   Total_IDR   = 0.0; Total_IDC           = 0.0;
  Total_Custom_ObjFunc = 0.0;
  
  /*--- Read farfield conditions from config ---*/
  
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Mach_Inf        = config->GetMach();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();
  Tke_Inf         = config->GetTke_FreeStreamND();
  
  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  
  switch(direct_diff) {
    case NO_DERIVATIVE:
      break;
    case D_DENSITY:
      SU2_TYPE::SetDerivative(Density_Inf, 1.0);
      break;
    case D_PRESSURE:
      SU2_TYPE::SetDerivative(Pressure_Inf, 1.0);
      break;
    case D_TEMPERATURE:
      SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
      break;
    case D_VISCOSITY:
      SU2_TYPE::SetDerivative(Viscosity_Inf, 1.0);
      break;
    case D_MACH: case D_AOA:
    case D_SIDESLIP: case D_REYNOLDS:
    case D_TURB2LAM: case D_DESIGN:
      /*--- Already done in postprocessing of config ---*/
      break;
    default:
      break;
  }

  /*--- Initialize fan face pressure, fan face mach number, and mass flow rate ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inflow_MassFlow[iMarker]     = 0.0;
    Inflow_Mach[iMarker]         = Mach_Inf;
    Inflow_Pressure[iMarker]     = Pressure_Inf;
    Inflow_Area[iMarker]         = 0.0;
    
    Exhaust_MassFlow[iMarker]    = 0.0;
    Exhaust_Temperature[iMarker] = Temperature_Inf;
    Exhaust_Pressure[iMarker]    = Pressure_Inf;
    Exhaust_Area[iMarker]        = 0.0;

  }
  /*--- Initializate quantities for SlidingMesh Interface ---*/
  
  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++){

    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;
    
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];

      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];

        SlidingStateNodes[iMarker][iPoint] = 0;
        for (iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }

    }
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/

  counter_local = 0;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    Density = node[iPoint]->GetSolution(0);

    Velocity2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);

    StaticEnergy= node[iPoint]->GetSolution(nDim+1)/Density - 0.5*Velocity2;

    FluidModel->SetTDState_rhoe(Density, StaticEnergy);
    Pressure= FluidModel->GetPressure();
    Temperature= FluidModel->GetTemperature();

    /*--- Use the values at the infinity ---*/

    if ((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
      Solution[0] = Density_Inf;
      for (iDim = 0; iDim < nDim; iDim++)
        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
      Solution[nDim+1] = Energy_Inf*Density_Inf;
      node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetSolution_Old(Solution);
      counter_local++;
    }

  }

  /*--- Warning message about non-physical points ---*/

  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  }

  /*--- Define solver parameters needed for execution of destructor ---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) space_centered = true;
  else space_centered = false;

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;

  /*--- Perform the MPI communication of the solution ---*/

  Set_MPI_Solution(geometry, config);
  
}

CNSSolver::~CNSSolver(void) {

  unsigned short iMarker, iDim;

  if (CD_Visc != NULL)          delete [] CD_Visc;
  if (CL_Visc != NULL)          delete [] CL_Visc;
  if (CSF_Visc != NULL)         delete [] CSF_Visc;
  if (CFx_Visc != NULL)         delete [] CFx_Visc;
  if (CFy_Visc != NULL)         delete [] CFy_Visc;
  if (CFz_Visc != NULL)         delete [] CFz_Visc;
  if (CMx_Visc != NULL)         delete [] CMx_Visc;
  if (CMy_Visc != NULL)         delete [] CMy_Visc;
  if (CMz_Visc != NULL)         delete [] CMz_Visc;
  if (CoPx_Visc != NULL)        delete [] CoPx_Visc;
  if (CoPy_Visc != NULL)        delete [] CoPy_Visc;
  if (CoPz_Visc != NULL)        delete [] CoPz_Visc;
  if (CEff_Visc != NULL)        delete [] CEff_Visc;
  if (CMerit_Visc != NULL)      delete [] CMerit_Visc;
  if (CT_Visc != NULL)          delete [] CT_Visc;
  if (CQ_Visc != NULL)          delete [] CQ_Visc;
  if (HF_Visc != NULL)          delete [] HF_Visc;
  if (MaxHF_Visc != NULL)       delete [] MaxHF_Visc;
  if (ForceViscous != NULL)     delete [] ForceViscous;
  if (MomentViscous != NULL)    delete [] MomentViscous;
  
  if (Surface_CL_Visc != NULL)      delete [] Surface_CL_Visc;
  if (Surface_CD_Visc != NULL)      delete [] Surface_CD_Visc;
  if (Surface_CSF_Visc != NULL)     delete [] Surface_CSF_Visc;
  if (Surface_CEff_Visc != NULL)    delete [] Surface_CEff_Visc;
  if (Surface_CFx_Visc != NULL)     delete [] Surface_CFx_Visc;
  if (Surface_CFy_Visc != NULL)     delete [] Surface_CFy_Visc;
  if (Surface_CFz_Visc != NULL)     delete [] Surface_CFz_Visc;
  if (Surface_CMx_Visc != NULL)     delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != NULL)     delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != NULL)     delete [] Surface_CMz_Visc;
  if (Surface_HF_Visc != NULL)      delete [] Surface_HF_Visc;
  if (Surface_MaxHF_Visc != NULL)   delete [] Surface_MaxHF_Visc;
  
  if (CSkinFriction != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        delete [] CSkinFriction[iMarker][iDim];
      }
      delete [] CSkinFriction[iMarker];
    }
    delete [] CSkinFriction;
  }
  
}

void CNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint, ErrorCounter = 0;
  su2double StrainMag = 0.0, Omega = 0.0, *Vorticity;
    
  unsigned long ExtIter     = config->GetExtIter();
  bool cont_adjoint         = config->GetContinuous_Adjoint();
  bool disc_adjoint         = config->GetDiscrete_Adjoint();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool center               = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst           = center && config->GetKind_Centered_Flow() == JST;
  bool limiter_flow         = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_turb         = ((config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_adjflow      = (cont_adjoint && (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool fixed_cl             = config->GetFixed_CL_Mode();
  bool engine               = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk        = ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0));
  bool nearfield            = (config->GetnMarker_NearFieldBound() != 0);
  bool interface            = (config->GetnMarker_InterfaceBound() != 0);
  bool van_albada           = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
  unsigned short kind_row_dissipation = config->GetKind_RoeLowDiss();
  bool roe_low_dissipation  = (kind_row_dissipation != NO_ROELOWDISS) && (config->GetKind_Upwind_Flow() == ROE);

  /*--- Update the angle of attack at the far-field for fixed CL calculations. ---*/
  
  if (fixed_cl) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

  /*--- Compute the engine properties ---*/

  if (engine) { GetPower_Properties(geometry, config, iMesh, Output); }

  /*--- Compute the actuator disk properties and distortion levels ---*/

  if (actuator_disk) {
    Set_MPI_ActDisk(solver_container, geometry, config);
    SetActDisk_BCThrust(geometry, solver_container, config, iMesh, Output);
  }

  /*--- Compute Interface MPI ---*/

  if (interface) { Set_MPI_Interface(geometry, config); }

  /*--- Compute NearField MPI ---*/

  if (nearfield) { Set_MPI_Nearfield(geometry, config); }
 
  /*--- Artificial dissipation ---*/

  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
  /*--- Roe Low Dissipation Sensor ---*/
  
  if (roe_low_dissipation){
    SetRoe_Dissipation(geometry, config);
    if (kind_row_dissipation == FD_DUCROS || kind_row_dissipation == NTS_DUCROS){
      SetUpwind_Ducros_Sensor(geometry, config);
    }
  }
  
  /*--- Compute gradient of the primitive variables ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
  }

  /*--- Compute the limiter in case we need it in the turbulence model
   or to limit the viscous terms (check this logic with JST and 2nd order turbulence model) ---*/

  if ((iMesh == MESH_0) && (limiter_flow || limiter_turb || limiter_adjflow)
      && !Output && !van_albada) { SetPrimitive_Limiter(geometry, config); }
  
  /*--- Evaluate the vorticity and strain rate magnitude ---*/
  
  StrainMag_Max = 0.0; Omega_Max = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    solver_container[FLOW_SOL]->node[iPoint]->SetVorticity();
    solver_container[FLOW_SOL]->node[iPoint]->SetStrainMag();
    
    StrainMag = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
    Vorticity = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity();
    Omega = sqrt(Vorticity[0]*Vorticity[0]+ Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);
    
    StrainMag_Max = max(StrainMag_Max, StrainMag);
    Omega_Max = max(Omega_Max, Omega);
    
  }

  /*--- Initialize the Jacobian matrices ---*/
  
  if (implicit && !config->GetDiscrete_Adjoint()) Jacobian.SetValZero();

  /*--- Error message ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    su2double MyOmega_Max = Omega_Max; Omega_Max = 0.0;
    su2double MyStrainMag_Max = StrainMag_Max; StrainMag_Max = 0.0;
    
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyStrainMag_Max, &StrainMag_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyOmega_Max, &Omega_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    
    if (iMesh == MESH_0) {
      config->SetNonphysical_Points(ErrorCounter);
      solver_container[FLOW_SOL]->SetStrainMag_Max(StrainMag_Max);
      solver_container[FLOW_SOL]->SetOmega_Max(Omega_Max);
    }
    
  }
  
}

unsigned long CNSSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double eddy_visc = 0.0, turb_ke = 0.0, DES_LengthScale = 0.0;
  unsigned short turb_model = config->GetKind_Turb_Model();
  bool RightSol = true;
  
  bool tkeNeeded            = (turb_model == SST);

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Retrieve the value of the kinetic energy (if need it) ---*/
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->node[iPoint]->GetmuT();
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
      
      if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
        DES_LengthScale = solver_container[TURB_SOL]->node[iPoint]->GetDES_LengthScale();
      }
    }
    
    /*--- Initialize the non-physical points vector ---*/
    
    node[iPoint]->SetNon_Physical(false);
    
    /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/
    
    RightSol = node[iPoint]->SetPrimVar(eddy_visc, turb_ke, FluidModel);
    node[iPoint]->SetSecondaryVar(FluidModel);

    if (!RightSol) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }
        
    /*--- Set the DES length scale ---*/
    
    node[iPoint]->SetDES_LengthScale(DES_LengthScale);
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  return ErrorCounter;
}

void CNSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
  
  su2double *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time, Local_Delta_Time_Visc,
  Global_Delta_Time = 1E6, Mean_LaminarVisc = 0.0, Mean_EddyVisc = 0.0, Mean_Density = 0.0, Lambda_1, Lambda_2, K_v = 0.25, Global_Delta_UnstTimeND;
  unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
  unsigned short iDim, iMarker;
  su2double ProjVel, ProjVel_i, ProjVel_j;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetMax_Lambda_Inv(0.0);
    node[iPoint]->SetMax_Lambda_Visc(0.0);
  }
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    
    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j) ;
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed ;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);
    
    /*--- Viscous contribution ---*/
    
    Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
    Mean_EddyVisc    = 0.5*(node[iPoint]->GetEddyViscosity() + node[jPoint]->GetEddyViscosity());
    Mean_Density     = 0.5*(node[iPoint]->GetSolution(0) + node[jPoint]->GetSolution(0));
    
    Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
    //TODO (REAL_GAS) removing Gamma it cannot work with FLUIDPROP
    Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
    Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
    
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Visc(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      
      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
      
      /*--- Viscous contribution ---*/
      
      Mean_LaminarVisc = node[iPoint]->GetLaminarViscosity();
      Mean_EddyVisc    = node[iPoint]->GetEddyViscosity();
      Mean_Density     = node[iPoint]->GetSolution(0);
      
      Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
      Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
      Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
      
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
      
    }
  }
  
  /*--- Each element uses their own speed, steady state simulation ---*/
  
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();
      Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
    
  }
  
  
  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;
    
    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
#endif
  }
  
  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  }
  
  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
    
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = rbuf_time;
#endif
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }
  
  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  
}

void CNSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iPoint, jPoint, iEdge;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points, coordinates and normal vector in edge ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive and secondary variables ---*/
    
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
    numerics->SetSecondary(node[iPoint]->GetSecondary(), node[jPoint]->GetSecondary());
    
    /*--- Gradient and limiters ---*/
    
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());
    
    /*--- Turbulent kinetic energy ---*/
    
    if (config->GetKind_Turb_Model() == SST)
      numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                     solver_container[TURB_SOL]->node[jPoint]->GetSolution(0));
    
    /*--- Compute and update residual ---*/
    
    numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.SubtractBlock(iPoint, Res_Visc);
    LinSysRes.AddBlock(jPoint, Res_Visc);
    
    /*--- Implicit part ---*/
    
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    }
    
  }
  
}

void CNSSolver::Friction_Forces(CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint, iPointNormal;
  unsigned short Boundary, Monitoring, iMarker, iMarker_Monitoring, iDim, jDim;
  su2double Viscosity = 0.0, div_vel, *Normal, MomentDist[3] = {0.0, 0.0, 0.0}, WallDist[3] = {0.0, 0.0, 0.0},
  *Coord, *Coord_Normal, Area, WallShearStress, TauNormal, factor, RefTemp, RefVel2,
  RefDensity, GradTemperature, Density = 0.0, WallDistMod, FrictionVel,
  Mach2Vel, Mach_Motion, UnitNormal[3] = {0.0, 0.0, 0.0}, TauElem[3] = {0.0, 0.0, 0.0}, TauTangent[3] = {0.0, 0.0, 0.0},
  Tau[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}}, Force[3] = {0.0, 0.0, 0.0}, Cp, thermal_conductivity, MaxNorm = 8.0,
  Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}}, Grad_Temp[3] = {0.0, 0.0, 0.0},
  delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  su2double MomentX_Force[3] = {0.0,0.0,0.0}, MomentY_Force[3] = {0.0,0.0,0.0}, MomentZ_Force[3] = {0.0,0.0,0.0};
  su2double AxiFactor;

#ifdef HAVE_MPI
  su2double MyAllBound_CD_Visc, MyAllBound_CL_Visc, MyAllBound_CSF_Visc, MyAllBound_CMx_Visc, MyAllBound_CMy_Visc, MyAllBound_CMz_Visc, MyAllBound_CoPx_Visc, MyAllBound_CoPy_Visc, MyAllBound_CoPz_Visc, MyAllBound_CFx_Visc, MyAllBound_CFy_Visc, MyAllBound_CFz_Visc, MyAllBound_CT_Visc, MyAllBound_CQ_Visc, MyAllBound_HF_Visc, MyAllBound_MaxHF_Visc, *MySurface_CL_Visc = NULL, *MySurface_CD_Visc = NULL, *MySurface_CSF_Visc = NULL, *MySurface_CEff_Visc = NULL, *MySurface_CFx_Visc = NULL, *MySurface_CFy_Visc = NULL, *MySurface_CFz_Visc = NULL, *MySurface_CMx_Visc = NULL, *MySurface_CMy_Visc = NULL, *MySurface_CMz_Visc = NULL, *MySurface_HF_Visc = NULL, *MySurface_MaxHF_Visc;
#endif
  
  string Marker_Tag, Monitoring_Tag;
  
  su2double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta            = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea    = config->GetRefArea();
  su2double RefLength = config->GetRefLength();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  su2double *Origin = NULL;
  
  if (config->GetnMarker_Monitoring() != 0) { Origin = config->GetRefOriginMoment(0); }
  
  bool grid_movement        = config->GetGrid_Movement();
  su2double Prandtl_Lam     = config->GetPrandtl_Lam();
  bool QCR                  = config->GetQCR();
  bool axisymmetric         = config->GetAxisymmetric();

  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/
  
  RefTemp    = Temperature_Inf;
  RefDensity = Density_Inf;
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  } else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);
  
  /*--- Variables initialization ---*/
  
  AllBound_CD_Visc = 0.0;    AllBound_CL_Visc = 0.0;       AllBound_CSF_Visc = 0.0;
  AllBound_CFx_Visc = 0.0;      AllBound_CFy_Visc = 0.0;         AllBound_CFz_Visc = 0.0;
  AllBound_CMx_Visc = 0.0;      AllBound_CMy_Visc = 0.0;         AllBound_CMz_Visc = 0.0;
  AllBound_CoPx_Visc = 0.0;      AllBound_CoPy_Visc = 0.0;         AllBound_CoPz_Visc = 0.0;
  AllBound_CT_Visc = 0.0;       AllBound_CQ_Visc = 0.0;          AllBound_CMerit_Visc = 0.0;
  AllBound_HF_Visc = 0.0; AllBound_MaxHF_Visc = 0.0; AllBound_CEff_Visc = 0.0;
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Visc[iMarker_Monitoring]      = 0.0; Surface_CD_Visc[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Visc[iMarker_Monitoring] = 0.0; Surface_CEff_Visc[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring]        = 0.0; Surface_CFy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring]        = 0.0; Surface_CMx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]        = 0.0; Surface_CMz_Visc[iMarker_Monitoring]        = 0.0;
    Surface_HF_Visc[iMarker_Monitoring]              = 0.0; Surface_MaxHF_Visc[iMarker_Monitoring]           = 0.0;
  }
  
  /*--- Loop over the Navier-Stokes markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Boundary = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    /*--- Obtain the origin for the moment computation for a particular marker ---*/
    
    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }
    
    if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {
      
      /*--- Forces initialization at each Marker ---*/
      
      CD_Visc[iMarker] = 0.0; CL_Visc[iMarker] = 0.0;       CSF_Visc[iMarker] = 0.0;
      CFx_Visc[iMarker] = 0.0;   CFy_Visc[iMarker] = 0.0;         CFz_Visc[iMarker] = 0.0;
      CMx_Visc[iMarker] = 0.0;   CMy_Visc[iMarker] = 0.0;         CMz_Visc[iMarker] = 0.0;
      CoPx_Visc[iMarker] = 0.0;  CoPy_Visc[iMarker] = 0.0;       CoPz_Visc[iMarker] = 0.0;
      CT_Visc[iMarker] = 0.0;    CQ_Visc[iMarker] = 0.0;          CMerit_Visc[iMarker] = 0.0;
      HF_Visc[iMarker] = 0.0;  MaxHF_Visc[iMarker] = 0.0; CEff_Visc[iMarker] = 0.0;
      
      for (iDim = 0; iDim < nDim; iDim++) ForceViscous[iDim] = 0.0;
      MomentViscous[0] = 0.0; MomentViscous[1] = 0.0; MomentViscous[2] = 0.0;
      MomentX_Force[0] = 0.0; MomentX_Force[1] = 0.0; MomentX_Force[2] = 0.0;
      MomentY_Force[0] = 0.0; MomentY_Force[1] = 0.0; MomentY_Force[2] = 0.0;
      MomentZ_Force[0] = 0.0; MomentZ_Force[1] = 0.0; MomentZ_Force[2] = 0.0;

      /*--- Loop over the vertices to compute the forces ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
        
        Coord = geometry->node[iPoint]->GetCoord();
        Coord_Normal = geometry->node[iPointNormal]->GetCoord();
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
          }
          Grad_Temp[iDim] = node[iPoint]->GetGradient_Primitive(0, iDim);
        }
        
        Viscosity = node[iPoint]->GetLaminarViscosity();
        Density = node[iPoint]->GetDensity();
        
        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
        

        for (iDim = 0; iDim < nDim; iDim++) {
          UnitNormal[iDim] = Normal[iDim]/Area;
        }
        
        /*--- Evaluate Tau ---*/
        
        div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim]) - TWO3*Viscosity*div_vel*delta[iDim][jDim];
          }
        }
        
        /*--- If necessary evaluate the QCR contribution to Tau ---*/
        
        if (QCR){
            su2double den_aux, c_cr1=0.3, O_ik, O_jk;
            unsigned short kDim;
            
            /*--- Denominator Antisymmetric normalized rotation tensor ---*/
            
            den_aux = 0.0;
            for (iDim = 0 ; iDim < nDim; iDim++)
                for (jDim = 0 ; jDim < nDim; jDim++)
                    den_aux += Grad_Vel[iDim][jDim] * Grad_Vel[iDim][jDim];
            den_aux = sqrt(max(den_aux,1E-10));
            
            /*--- Adding the QCR contribution ---*/
            
            for (iDim = 0 ; iDim < nDim; iDim++){
                for (jDim = 0 ; jDim < nDim; jDim++){
                    for (kDim = 0 ; kDim < nDim; kDim++){
                        O_ik = (Grad_Vel[iDim][kDim] - Grad_Vel[kDim][iDim])/ den_aux;
                        O_jk = (Grad_Vel[jDim][kDim] - Grad_Vel[kDim][jDim])/ den_aux;
                        Tau[iDim][jDim] -= c_cr1 * (O_ik * Tau[jDim][kDim] + O_jk * Tau[iDim][kDim]);
                    }
                }
            }
        
        }
        
        /*--- Project Tau in each surface element ---*/
        
        for (iDim = 0; iDim < nDim; iDim++) {
          TauElem[iDim] = 0.0;
          for (jDim = 0; jDim < nDim; jDim++) {
            TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
          }
        }
        
        /*--- Compute wall shear stress (using the stress tensor). Compute wall skin friction coefficient, and heat flux on the wall ---*/
        
        TauNormal = 0.0; for (iDim = 0; iDim < nDim; iDim++) TauNormal += TauElem[iDim] * UnitNormal[iDim];
        
        WallShearStress = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
          CSkinFriction[iMarker][iDim][iVertex] = TauTangent[iDim] / (0.5*RefDensity*RefVel2);
          WallShearStress += TauTangent[iDim] * TauTangent[iDim];
        }
        WallShearStress = sqrt(WallShearStress);
        
        for (iDim = 0; iDim < nDim; iDim++) WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
        WallDistMod = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim]; WallDistMod = sqrt(WallDistMod);
        
        /*--- Compute y+ and non-dimensional velocity ---*/
        
        FrictionVel = sqrt(fabs(WallShearStress)/Density);
        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);
        
        /*--- Compute total and maximum heat flux on the wall ---*/

        GradTemperature = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          GradTemperature -= Grad_Temp[iDim]*UnitNormal[iDim];

        Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
        thermal_conductivity = Cp * Viscosity/Prandtl_Lam;
        HeatFlux[iMarker][iVertex] = -thermal_conductivity*GradTemperature;
        HF_Visc[iMarker] += HeatFlux[iMarker][iVertex]*Area;
        MaxHF_Visc[iMarker] += pow(HeatFlux[iMarker][iVertex], MaxNorm);

        
        /*--- Note that y+, and heat are computed at the
         halo cells (for visualization purposes), but not the forces ---*/
        
        if ((geometry->node[iPoint]->GetDomain()) && (Monitoring == YES)) {
          
          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Force computation ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = TauElem[iDim] * Area * factor * AxiFactor;
            ForceViscous[iDim] += Force[iDim];
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }
          
          /*--- Moment with respect to the reference axis ---*/
          
          if (iDim == 3) {
            MomentViscous[0] += (Force[2]*MomentDist[1] - Force[1]*MomentDist[2])/RefLength;
            MomentX_Force[1] += (-Force[1]*Coord[2]);
            MomentX_Force[2] += (Force[2]*Coord[1]);

            MomentViscous[1] += (Force[0]*MomentDist[2] - Force[2]*MomentDist[0])/RefLength;
            MomentY_Force[2] += (-Force[2]*Coord[0]);
            MomentY_Force[0] += (Force[0]*Coord[2]);
          }
          MomentViscous[2] += (Force[1]*MomentDist[0] - Force[0]*MomentDist[1])/RefLength;
          MomentZ_Force[0] += (-Force[0]*Coord[1]);
          MomentZ_Force[1] += (Force[1]*Coord[0]);
          
        }
        
      }
      
      /*--- Project forces and store the non-dimensional coefficients ---*/
      
      if (Monitoring == YES) {
        if (nDim == 2) {
          CD_Visc[iMarker]          =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          CL_Visc[iMarker]          = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          CEff_Visc[iMarker]        = CL_Visc[iMarker] / (CD_Visc[iMarker]+EPS);
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CMz_Visc[iMarker]         = MomentViscous[2];
          CoPx_Visc[iMarker]        = MomentZ_Force[1];
          CoPy_Visc[iMarker]        = -MomentZ_Force[0];
          CT_Visc[iMarker]          = -CFx_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker] / (CQ_Visc[iMarker]+EPS);
          MaxHF_Visc[iMarker]       = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }
        if (nDim == 3) {
          CD_Visc[iMarker]          =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
          CL_Visc[iMarker]          = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          CSF_Visc[iMarker]         = -ForceViscous[0]*sin(Beta)*cos(Alpha) + ForceViscous[1]*cos(Beta) - ForceViscous[2]*sin(Beta)*sin(Alpha);
          CEff_Visc[iMarker]        = CL_Visc[iMarker]/(CD_Visc[iMarker] + EPS);
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CFz_Visc[iMarker]         = ForceViscous[2];
          CMx_Visc[iMarker]         = MomentViscous[0];
          CMy_Visc[iMarker]         = MomentViscous[1];
          CMz_Visc[iMarker]         = MomentViscous[2];
          CoPx_Visc[iMarker]        =  -MomentY_Force[0];
          CoPz_Visc[iMarker]        = MomentY_Force[2];
          CT_Visc[iMarker]          = -CFz_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker] / (CQ_Visc[iMarker] + EPS);
          MaxHF_Visc[iMarker]       = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }
        
        AllBound_CD_Visc          += CD_Visc[iMarker];
        AllBound_CL_Visc          += CL_Visc[iMarker];
        AllBound_CSF_Visc         += CSF_Visc[iMarker];
        AllBound_CFx_Visc         += CFx_Visc[iMarker];
        AllBound_CFy_Visc         += CFy_Visc[iMarker];
        AllBound_CFz_Visc         += CFz_Visc[iMarker];
        AllBound_CMx_Visc         += CMx_Visc[iMarker];
        AllBound_CMy_Visc         += CMy_Visc[iMarker];
        AllBound_CMz_Visc         += CMz_Visc[iMarker];
        AllBound_CoPx_Visc        += CoPx_Visc[iMarker];
        AllBound_CoPy_Visc        += CoPy_Visc[iMarker];
        AllBound_CoPz_Visc        += CoPz_Visc[iMarker];
        AllBound_CT_Visc          += CT_Visc[iMarker];
        AllBound_CQ_Visc          += CQ_Visc[iMarker];
        AllBound_HF_Visc          += HF_Visc[iMarker];
        AllBound_MaxHF_Visc       += pow(MaxHF_Visc[iMarker], MaxNorm);
        
        /*--- Compute the coefficients per surface ---*/
        
        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Visc[iMarker_Monitoring]      += CL_Visc[iMarker];
            Surface_CD_Visc[iMarker_Monitoring]      += CD_Visc[iMarker];
            Surface_CSF_Visc[iMarker_Monitoring] += CSF_Visc[iMarker];
            Surface_CEff_Visc[iMarker_Monitoring]       += CEff_Visc[iMarker];
            Surface_CFx_Visc[iMarker_Monitoring]        += CFx_Visc[iMarker];
            Surface_CFy_Visc[iMarker_Monitoring]        += CFy_Visc[iMarker];
            Surface_CFz_Visc[iMarker_Monitoring]        += CFz_Visc[iMarker];
            Surface_CMx_Visc[iMarker_Monitoring]        += CMx_Visc[iMarker];
            Surface_CMy_Visc[iMarker_Monitoring]        += CMy_Visc[iMarker];
            Surface_CMz_Visc[iMarker_Monitoring]        += CMz_Visc[iMarker];
            Surface_HF_Visc[iMarker_Monitoring]         += HF_Visc[iMarker];
            Surface_MaxHF_Visc[iMarker_Monitoring]      += pow(MaxHF_Visc[iMarker],MaxNorm);
          }
        }
        
      }
      
    }
  }
  
  /*--- Update some global coeffients ---*/
  
  AllBound_CEff_Visc = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  AllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);
  
  
#ifdef HAVE_MPI
  
  /*--- Add AllBound information using all the nodes ---*/
  
  MyAllBound_CD_Visc        = AllBound_CD_Visc;                      AllBound_CD_Visc = 0.0;
  MyAllBound_CL_Visc        = AllBound_CL_Visc;                      AllBound_CL_Visc = 0.0;
  MyAllBound_CSF_Visc   = AllBound_CSF_Visc;                 AllBound_CSF_Visc = 0.0;
  AllBound_CEff_Visc = 0.0;
  MyAllBound_CMx_Visc          = AllBound_CMx_Visc;                        AllBound_CMx_Visc = 0.0;
  MyAllBound_CMy_Visc          = AllBound_CMy_Visc;                        AllBound_CMy_Visc = 0.0;
  MyAllBound_CMz_Visc          = AllBound_CMz_Visc;                        AllBound_CMz_Visc = 0.0;
  MyAllBound_CoPx_Visc          = AllBound_CoPx_Visc;                        AllBound_CoPx_Visc = 0.0;
  MyAllBound_CoPy_Visc          = AllBound_CoPy_Visc;                        AllBound_CoPy_Visc = 0.0;
  MyAllBound_CoPz_Visc          = AllBound_CoPz_Visc;                        AllBound_CoPz_Visc = 0.0;
  MyAllBound_CFx_Visc          = AllBound_CFx_Visc;                        AllBound_CFx_Visc = 0.0;
  MyAllBound_CFy_Visc          = AllBound_CFy_Visc;                        AllBound_CFy_Visc = 0.0;
  MyAllBound_CFz_Visc          = AllBound_CFz_Visc;                        AllBound_CFz_Visc = 0.0;
  MyAllBound_CT_Visc           = AllBound_CT_Visc;                         AllBound_CT_Visc = 0.0;
  MyAllBound_CQ_Visc           = AllBound_CQ_Visc;                         AllBound_CQ_Visc = 0.0;
  AllBound_CMerit_Visc = 0.0;
  MyAllBound_HF_Visc     = AllBound_HF_Visc;                   AllBound_HF_Visc = 0.0;
  MyAllBound_MaxHF_Visc  = pow(AllBound_MaxHF_Visc, MaxNorm);  AllBound_MaxHF_Visc = 0.0;
  
  SU2_MPI::Allreduce(&MyAllBound_CD_Visc, &AllBound_CD_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CL_Visc, &AllBound_CL_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CSF_Visc, &AllBound_CSF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Visc = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CMx_Visc, &AllBound_CMx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMy_Visc, &AllBound_CMy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMz_Visc, &AllBound_CMz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFx_Visc, &AllBound_CFx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFy_Visc, &AllBound_CFy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFz_Visc, &AllBound_CFz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CoPx_Visc, &AllBound_CoPx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CoPy_Visc, &AllBound_CoPy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CoPz_Visc, &AllBound_CoPz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CT_Visc, &AllBound_CT_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CQ_Visc, &AllBound_CQ_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  SU2_MPI::Allreduce(&MyAllBound_HF_Visc, &AllBound_HF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_MaxHF_Visc, &AllBound_MaxHF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);
  
  /*--- Add the forces on the surfaces using all the nodes ---*/
  
  MySurface_CL_Visc         = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Visc         = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Visc       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_HF_Visc         = new su2double[config->GetnMarker_Monitoring()];
  MySurface_MaxHF_Visc      = new su2double[config->GetnMarker_Monitoring()];
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    
    MySurface_CL_Visc[iMarker_Monitoring]      = Surface_CL_Visc[iMarker_Monitoring];
    MySurface_CD_Visc[iMarker_Monitoring]      = Surface_CD_Visc[iMarker_Monitoring];
    MySurface_CSF_Visc[iMarker_Monitoring] = Surface_CSF_Visc[iMarker_Monitoring];
    MySurface_CEff_Visc[iMarker_Monitoring]       = Surface_CEff_Visc[iMarker_Monitoring];
    MySurface_CFx_Visc[iMarker_Monitoring]        = Surface_CFx_Visc[iMarker_Monitoring];
    MySurface_CFy_Visc[iMarker_Monitoring]        = Surface_CFy_Visc[iMarker_Monitoring];
    MySurface_CFz_Visc[iMarker_Monitoring]        = Surface_CFz_Visc[iMarker_Monitoring];
    MySurface_CMx_Visc[iMarker_Monitoring]        = Surface_CMx_Visc[iMarker_Monitoring];
    MySurface_CMy_Visc[iMarker_Monitoring]        = Surface_CMy_Visc[iMarker_Monitoring];
    MySurface_CMz_Visc[iMarker_Monitoring]        = Surface_CMz_Visc[iMarker_Monitoring];
    MySurface_HF_Visc[iMarker_Monitoring]         = Surface_HF_Visc[iMarker_Monitoring];
    MySurface_MaxHF_Visc[iMarker_Monitoring]      = Surface_MaxHF_Visc[iMarker_Monitoring];
    
    Surface_CL_Visc[iMarker_Monitoring]         = 0.0;
    Surface_CD_Visc[iMarker_Monitoring]         = 0.0;
    Surface_CSF_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CEff_Visc[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Visc[iMarker_Monitoring]        = 0.0;
    Surface_HF_Visc[iMarker_Monitoring]         = 0.0;
    Surface_MaxHF_Visc[iMarker_Monitoring]      = 0.0;
  }
  
  SU2_MPI::Allreduce(MySurface_CL_Visc, Surface_CL_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CD_Visc, Surface_CD_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CSF_Visc, Surface_CSF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
    Surface_CEff_Visc[iMarker_Monitoring] = Surface_CL_Visc[iMarker_Monitoring] / (Surface_CD_Visc[iMarker_Monitoring] + EPS);
  SU2_MPI::Allreduce(MySurface_CFx_Visc, Surface_CFx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFy_Visc, Surface_CFy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFz_Visc, Surface_CFz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMx_Visc, Surface_CMx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMy_Visc, Surface_CMy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMz_Visc, Surface_CMz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_HF_Visc, Surface_HF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_MaxHF_Visc, Surface_MaxHF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  delete [] MySurface_CL_Visc; delete [] MySurface_CD_Visc; delete [] MySurface_CSF_Visc;
  delete [] MySurface_CEff_Visc;  delete [] MySurface_CFx_Visc;   delete [] MySurface_CFy_Visc;
  delete [] MySurface_CFz_Visc;   delete [] MySurface_CMx_Visc;   delete [] MySurface_CMy_Visc;
  delete [] MySurface_CMz_Visc;   delete [] MySurface_HF_Visc; delete [] MySurface_MaxHF_Visc;
  
#endif
  
  /*--- Update the total coefficients (note that all the nodes have the same value)---*/
  
  Total_CD          += AllBound_CD_Visc;
  Total_CL          += AllBound_CL_Visc;
  Total_CSF         += AllBound_CSF_Visc;
  Total_CEff        = Total_CL / (Total_CD + EPS);
  Total_CFx         += AllBound_CFx_Visc;
  Total_CFy         += AllBound_CFy_Visc;
  Total_CFz         += AllBound_CFz_Visc;
  Total_CMx         += AllBound_CMx_Visc;
  Total_CMy         += AllBound_CMy_Visc;
  Total_CMz         += AllBound_CMz_Visc;
  Total_CoPx        += AllBound_CoPx_Visc;
  Total_CoPy        += AllBound_CoPy_Visc;
  Total_CoPz        += AllBound_CoPz_Visc;
  Total_CT          += AllBound_CT_Visc;
  Total_CQ          += AllBound_CQ_Visc;
  Total_CMerit      = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  Total_Heat        = AllBound_HF_Visc;
  Total_MaxHeat     = AllBound_MaxHF_Visc;
  
  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]      += Surface_CL_Visc[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]      += Surface_CD_Visc[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring] += Surface_CSF_Visc[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       = Surface_CL[iMarker_Monitoring] / (Surface_CD[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        += Surface_CFx_Visc[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        += Surface_CFy_Visc[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        += Surface_CFz_Visc[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        += Surface_CMx_Visc[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        += Surface_CMy_Visc[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        += Surface_CMz_Visc[iMarker_Monitoring];
  }
  
}

void CNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, jDim, iVar, jVar, Wall_Function;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  
  su2double Wall_HeatFlux, dist_ij, *Coord_i, *Coord_j, theta2;
  su2double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  su2double ProjGridVel, *GridVel, GridVel2, *Normal, Area, Pressure = 0.0;
  su2double total_viscosity, div_vel, Density, tau_vel[3] = {0.0, 0.0, 0.0}, UnitNormal[3] = {0.0, 0.0, 0.0};
  su2double laminar_viscosity = 0.0, eddy_viscosity = 0.0, Grad_Vel[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}},
  tau[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  su2double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Identify the boundary by string name ---*/
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Get the specified wall heat flux from config as well as the
        wall function treatment.---*/
  
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);
  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
  if(Wall_Function != NO_WALL_FUNCTION) {
    SU2_MPI::Error("Wall function treament not implemented yet", CURRENT_FUNCTION);
  }
  
  /*--- Loop over all of the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {

      /* If it is a customizable patch, retrieve the specified wall heat flux. */
      if (config->GetMarker_All_PyCustom(val_marker)) Wall_HeatFlux = geometry->GetCustomBoundaryHeatFlux(val_marker, iVertex);
      
      /*--- Compute dual-grid area and boundary normal ---*/
      
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }
      
      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there are moving walls (v = u_wall)---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/
      
      node[iPoint]->SetVelocity_Old(Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      
      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/
      
      Res_Visc[nDim+1] = Wall_HeatFlux * Area;
      
      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/
      
      if (grid_movement) {
        
        /*--- Get the grid velocity at the current boundary node ---*/
        
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Retrieve other primitive quantities and viscosities ---*/
        
        Density  = node[iPoint]->GetSolution(0);
        Pressure = node[iPoint]->GetPressure();
        laminar_viscosity = node[iPoint]->GetLaminarViscosity();
        eddy_viscosity    = node[iPoint]->GetEddyViscosity();
        total_viscosity   = laminar_viscosity + eddy_viscosity;
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
          }
        }
        
        /*--- Divergence of the velocity ---*/
        
        div_vel = 0.0; for (iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        /*--- Compute the viscous stress tensor ---*/
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim]+Grad_Vel[iDim][jDim] ) - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
          }
        }
        
        /*--- Dot product of the stress tensor with the grid velocity ---*/
        
        for (iDim = 0 ; iDim < nDim; iDim++) {
          tau_vel[iDim] = 0.0;
          for (jDim = 0 ; jDim < nDim; jDim++)
            tau_vel[iDim] += tau[iDim][jDim]*GridVel[jDim];
        }
        
        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/
        
        Res_Conv[nDim+1] = Pressure*ProjGridVel;
        for (iDim = 0 ; iDim < nDim; iDim++)
          Res_Visc[nDim+1] += tau_vel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Implicit Jacobian contributions due to moving walls ---*/
        
        if (implicit) {
          
          /*--- Jacobian contribution related to the pressure term ---*/
          
          GridVel2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            GridVel2 += GridVel[iDim]*GridVel[iDim];
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_i[nDim+1][0] = 0.5*(Gamma-1.0)*GridVel2*ProjGridVel;
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nDim+1][jDim+1] = -(Gamma-1.0)*GridVel[jDim]*ProjGridVel;
          Jacobian_i[nDim+1][nDim+1] = (Gamma-1.0)*ProjGridVel;
          
          /*--- Add the block to the Global Jacobian structure ---*/
          
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
          
          /*--- Now the Jacobian contribution related to the shear stress ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          
          /*--- Compute closest normal neighbor ---*/
          
          Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
          
          /*--- Get coordinates of i & nearest normal and compute distance ---*/
          
          Coord_i = geometry->node[iPoint]->GetCoord();
          Coord_j = geometry->node[Point_Normal]->GetCoord();
          
          dist_ij = 0;
          for (iDim = 0; iDim < nDim; iDim++)
            dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
          dist_ij = sqrt(dist_ij);
          
          theta2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            theta2 += UnitNormal[iDim]*UnitNormal[iDim];
          
          factor = total_viscosity*Area/(Density*dist_ij);
          
          if (nDim == 2) {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            
            etaz   = UnitNormal[0]*UnitNormal[1]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
          } else {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            thetaz = theta2 + UnitNormal[2]*UnitNormal[2]/3.0;
            
            etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            etax = UnitNormal[1]*UnitNormal[2]/3.0;
            etay = UnitNormal[0]*UnitNormal[2]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz   + GridVel[2]*etay;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay + GridVel[2]*etax;
            piz = GridVel[0]*etay   + GridVel[1]*etax   + GridVel[2]*thetaz;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]+piz*GridVel[2]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
            Jacobian_i[nDim+1][3] -= factor*piz;
          }
          
          /*--- Subtract the block from the Global Jacobian structure ---*/
          
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
          
        }
      }
      
      /*--- Convective contribution to the residual at the wall ---*/
      
      LinSysRes.AddBlock(iPoint, Res_Conv);
      
      /*--- Viscous contribution to the residual at the wall ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      
      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/
      
      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
      
    }
  }
}

void CNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iVar, jVar, iDim, jDim, Wall_Function;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  
  su2double *Normal, *Coord_i, *Coord_j, Area, dist_ij, theta2;
  su2double Twall, dTdn, dTdrho, thermal_conductivity;
  su2double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  su2double ProjGridVel, *GridVel, GridVel2, Pressure = 0.0, Density, Vel2;
  su2double total_viscosity, div_vel, tau_vel[3] = {0.0,0.0,0.0}, UnitNormal[3] = {0.0,0.0,0.0};
  su2double laminar_viscosity, eddy_viscosity, Grad_Vel[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}},
  tau[3][3] = {{0.0, 0.0, 0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}, delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
  su2double Prandtl_Lam  = config->GetPrandtl_Lam();
  su2double Prandtl_Turb = config->GetPrandtl_Turb();
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Identify the boundary ---*/
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Retrieve the specified wall temperature from config
        as well as the wall function treatment.---*/
  
  Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();
  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
  if(Wall_Function != NO_WALL_FUNCTION) {
    SU2_MPI::Error("Wall function treament not implemented yet", CURRENT_FUNCTION);
  }
  
  /*--- Loop over boundary points ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {

      /* If it is a customizable patch, retrieve the specified wall temperature. */
      if (config->GetMarker_All_PyCustom(val_marker)) Twall = geometry->GetCustomBoundaryTemperature(val_marker, iVertex);
      
      /*--- Compute dual-grid area and boundary normal ---*/
      
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Calculate useful quantities ---*/
      
      theta2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        theta2 += UnitNormal[iDim]*UnitNormal[iDim];
      
      /*--- Compute closest normal neighbor ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Get coordinates of i & nearest normal and compute distance ---*/
      
      Coord_i = geometry->node[iPoint]->GetCoord();
      Coord_j = geometry->node[Point_Normal]->GetCoord();
      dist_ij = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      dist_ij = sqrt(dist_ij);
      
      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there is grid motion (v = u_wall)---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      }
      else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }
      
      /*--- Set the residual, truncation error and velocity value on the boundary ---*/
      
      node[iPoint]->SetVelocity_Old(Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      
      /*--- Compute the normal gradient in temperature using Twall ---*/
      
      dTdn = -(node[Point_Normal]->GetPrimitive(0) - Twall)/dist_ij;
      
      /*--- Get transport coefficients ---*/
      
      laminar_viscosity    = node[iPoint]->GetLaminarViscosity();
      eddy_viscosity       = node[iPoint]->GetEddyViscosity();
      thermal_conductivity = Cp * ( laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb);
      
      // work in progress on real-gases...
      //thermal_conductivity = node[iPoint]->GetThermalConductivity();
      //Cp = node[iPoint]->GetSpecificHeatCp();
      //thermal_conductivity += Cp*eddy_viscosity/Prandtl_Turb;
      
      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/
      
      Res_Visc[nDim+1] = thermal_conductivity * dTdn * Area;
      
      /*--- Calculate Jacobian for implicit time stepping ---*/
      
      if (implicit) {
        
        for (iVar = 0; iVar < nVar; iVar ++)
          for (jVar = 0; jVar < nVar; jVar ++)
            Jacobian_i[iVar][jVar] = 0.0;
        
        /*--- Calculate useful quantities ---*/
        
        Density = node[iPoint]->GetPrimitive(nDim+2);
        Vel2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Vel2 += node[iPoint]->GetPrimitive(iDim+1) * node[iPoint]->GetPrimitive(iDim+1);
        dTdrho = 1.0/Density * ( -Twall + (Gamma-1.0)/Gas_Constant*(Vel2/2.0) );
        
        /*--- Enforce the no-slip boundary condition in a strong way ---*/
        
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
        
        /*--- Add contributions to the Jacobian from the weak enforcement of the energy equations ---*/
        
        Jacobian_i[nDim+1][0]      = -thermal_conductivity*theta2/dist_ij * dTdrho * Area;
        Jacobian_i[nDim+1][nDim+1] = -thermal_conductivity*theta2/dist_ij * (Gamma-1.0)/(Gas_Constant*Density) * Area;
        
        /*--- Subtract the block from the Global Jacobian structure ---*/
        
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/
      
      if (grid_movement) {
        
        /*--- Get the grid velocity at the current boundary node ---*/
        
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Retrieve other primitive quantities and viscosities ---*/
        
        Density  = node[iPoint]->GetSolution(0);
        Pressure = node[iPoint]->GetPressure();
        laminar_viscosity = node[iPoint]->GetLaminarViscosity();
        eddy_viscosity    = node[iPoint]->GetEddyViscosity();
        
        total_viscosity   = laminar_viscosity + eddy_viscosity;
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
          }
        }
        
        /*--- Divergence of the velocity ---*/
        
        div_vel = 0.0; for (iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        /*--- Compute the viscous stress tensor ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim] ) - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
          }
        
        /*--- Dot product of the stress tensor with the grid velocity ---*/
        
        for (iDim = 0 ; iDim < nDim; iDim++) {
          tau_vel[iDim] = 0.0;
          for (jDim = 0 ; jDim < nDim; jDim++)
            tau_vel[iDim] += tau[iDim][jDim]*GridVel[jDim];
        }
        
        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/
        
        Res_Conv[nDim+1] = Pressure*ProjGridVel;
        for (iDim = 0 ; iDim < nDim; iDim++)
          Res_Visc[nDim+1] += tau_vel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Implicit Jacobian contributions due to moving walls ---*/
        
        if (implicit) {
          
          /*--- Jacobian contribution related to the pressure term ---*/
          
          GridVel2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            GridVel2 += GridVel[iDim]*GridVel[iDim];
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          
          Jacobian_i[nDim+1][0] = 0.5*(Gamma-1.0)*GridVel2*ProjGridVel;
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nDim+1][jDim+1] = -(Gamma-1.0)*GridVel[jDim]*ProjGridVel;
          Jacobian_i[nDim+1][nDim+1] = (Gamma-1.0)*ProjGridVel;
          
          /*--- Add the block to the Global Jacobian structure ---*/
          
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
          
          /*--- Now the Jacobian contribution related to the shear stress ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          
          factor = total_viscosity*Area/(Density*dist_ij);
          
          if (nDim == 2) {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            
            etaz   = UnitNormal[0]*UnitNormal[1]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
          }
          else {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            thetaz = theta2 + UnitNormal[2]*UnitNormal[2]/3.0;
            
            etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            etax = UnitNormal[1]*UnitNormal[2]/3.0;
            etay = UnitNormal[0]*UnitNormal[2]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz   + GridVel[2]*etay;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay + GridVel[2]*etax;
            piz = GridVel[0]*etay   + GridVel[1]*etax   + GridVel[2]*thetaz;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]+piz*GridVel[2]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
            Jacobian_i[nDim+1][3] -= factor*piz;
          }
          
          /*--- Subtract the block from the Global Jacobian structure ---*/
          
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }
        
      }
      
      /*--- Convective contribution to the residual at the wall ---*/
      
      LinSysRes.AddBlock(iPoint, Res_Conv);
      
      /*--- Viscous contribution to the residual at the wall ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      
      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/
      
      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
      
    }
  }
}

void CNSSolver::SetRoe_Dissipation(CGeometry *geometry, CConfig *config){
  
  unsigned long iPoint;
  su2double wall_distance;
  
  unsigned short kind_roe_dissipation = config->GetKind_RoeLowDiss();
  
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    
    if (kind_roe_dissipation == FD || kind_roe_dissipation == FD_DUCROS){
      if (config->GetKind_HybridRANSLES() == NO_HYBRIDRANSLES){
        wall_distance = geometry->node[iPoint]->GetWall_Distance();
      } else {
        wall_distance = node[iPoint]->GetDES_LengthScale();
      }
      
      node[iPoint]->SetRoe_Dissipation_FD(wall_distance);
    }
    
    if (kind_roe_dissipation == NTS || kind_roe_dissipation == NTS_DUCROS){
      node[iPoint]->SetRoe_Dissipation_NTS();
    }
  }
}
