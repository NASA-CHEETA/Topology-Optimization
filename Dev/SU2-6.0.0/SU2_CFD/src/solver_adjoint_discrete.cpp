/*!
 * \file solver_adjoint_discrete.cpp
 * \brief Main subroutines for solving the discrete adjoint problem.
 * \author T. Albring
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

CDiscAdjSolver::CDiscAdjSolver(void) : CSolver () {

}

CDiscAdjSolver::CDiscAdjSolver(CGeometry *geometry, CConfig *config)  : CSolver() {

}

CDiscAdjSolver::CDiscAdjSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver, unsigned short Kind_Solver, unsigned short iMesh)  : CSolver() {

  unsigned short iVar, iMarker, iDim;
  unsigned long iVertex, iPoint;
  string text_line, mesh_filename;
  ifstream restart_file;
  string filename, AdjExt;

  bool fsi = config->GetFSI_Simulation();

  nVar = direct_solver->GetnVar();
  nDim = geometry->GetnDim();

  /*--- Initialize arrays to NULL ---*/

  CSensitivity = NULL;

  /*-- Store some information about direct solver ---*/
  this->KindDirect_Solver = Kind_Solver;
  this->direct_solver = direct_solver;


  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Allocate the node variables ---*/

  node = new CVariable*[nPoint];

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 1.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 1.0;

  /*--- Define some auxiliary vectors related to the geometry adjoint (nDim) ---*/
  Solution_Geometry = new su2double[nDim];     for (iDim = 0; iDim < nDim; iDim++) Solution_Geometry[iDim] = 1.0;

  /*--- Define some auxiliary vectors related to the residual for problems with a BGS strategy---*/

  if (fsi){

    Residual_BGS      = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]      = 1.0;
    Residual_Max_BGS  = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 1.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar] = 0;
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
    }

  }


  /*--- Define some structures for locating max residuals ---*/

  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution   = new su2double[nVar];

  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 1e-16;

  /*--- Sensitivity definition and coefficient in all the markers ---*/

  CSensitivity = new su2double* [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      CSensitivity[iMarker]        = new su2double [geometry->nVertex[iMarker]];
  }

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          CSensitivity[iMarker][iVertex] = 0.0;
      }
  }

  /*--- Initialize the discrete adjoint solution to zero everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CDiscAdjVariable(Solution, nDim, nVar, config);

}

CDiscAdjSolver::~CDiscAdjSolver(void) { 

  unsigned short iMarker;

  if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] CSensitivity[iMarker];
    }
    delete [] CSensitivity;
  }

}

void CDiscAdjSolver::SetRecording(CGeometry* geometry, CConfig *config){


  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  unsigned long iPoint;
  unsigned short iVar;

  /*--- Reset the solution to the initial (converged) solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->node[iPoint]->SetSolution(node[iPoint]->GetSolution_Direct());
  }

  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n()[iVar]);
      }
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n1()[iVar]);
      }
    }
  }

  /*--- Set the Jacobian to zero since this is not done inside the fluid iteration
   * when running the discrete adjoint solver. ---*/

  direct_solver->Jacobian.SetValZero();

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void CDiscAdjSolver::SetMesh_Recording(CGeometry** geometry, CVolumetricMovement *grid_movement, CConfig *config) {


//  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
//      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
//  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

//  unsigned long ExtIter = config->GetExtIter();

  unsigned long iPoint;
  unsigned short iDim;

  /*--- Reset the solution to the initial (converged) position ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iDim = 0; iDim < nDim; iDim++){
      geometry[MESH_0]->node[iPoint]->SetCoord(iDim,node[iPoint]->GetGeometry_Direct(iDim));
    }
  }

  /*--- After moving all nodes, update the dual mesh. Recompute the edges and
   dual mesh control volumes in the domain and on the boundaries. ---*/

  grid_movement->UpdateDualGrid(geometry[MESH_0], config);

  /*--- After updating the dual mesh, compute the grid velocities (only dynamic problems). ---*/
//  if (time_n_needed){
//    geometry[MESH_0]->SetGridVelocity(config, ExtIter);
//  }

  /*--- Update the multigrid structure after moving the finest grid,
   including computing the grid velocities on the coarser levels. ---*/

  grid_movement->UpdateMultiGrid(geometry, config);

//  if (time_n_needed){
//    for (iPoint = 0; iPoint < nPoint; iPoint++){
//      for (iVar = 0; iVar < nVar; iVar++){
//        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n()[iVar]);
//      }
//    }
//  }
//  if (time_n1_needed){
//    for (iPoint = 0; iPoint < nPoint; iPoint++){
//      for (iVar = 0; iVar < nVar; iVar++){
//        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n1()[iVar]);
//      }
//    }
//  }

}

void CDiscAdjSolver::RegisterSolution(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, nPoint = geometry->GetnPoint();

  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND,
  input = true;

  /*--- Register solution at all necessary time instances and other variables on the tape ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->node[iPoint]->RegisterSolution(input);
  }
  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      direct_solver->node[iPoint]->RegisterSolution_time_n();
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      direct_solver->node[iPoint]->RegisterSolution_time_n1();
    }
  }
}

void CDiscAdjSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset) {

  /*--- Register farfield values as input ---*/

  if((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS && !config->GetBoolTurbomachinery())) {

    su2double Velocity_Ref = config->GetVelocity_Ref();
    Alpha                  = config->GetAoA()*PI_NUMBER/180.0;
    Beta                   = config->GetAoS()*PI_NUMBER/180.0;
    Mach                   = config->GetMach();
    Pressure               = config->GetPressure_FreeStreamND();
    Temperature            = config->GetTemperature_FreeStreamND();

    su2double SoundSpeed = 0.0;
    
    if (nDim == 2) { SoundSpeed = config->GetVelocity_FreeStreamND()[0]*Velocity_Ref/(cos(Alpha)*Mach); }
    if (nDim == 3) { SoundSpeed = config->GetVelocity_FreeStreamND()[0]*Velocity_Ref/(cos(Alpha)*cos(Beta)*Mach); }

    if (!reset) {
      AD::RegisterInput(Mach);
      AD::RegisterInput(Alpha);
      AD::RegisterInput(Temperature);
      AD::RegisterInput(Pressure);
    }

    /*--- Recompute the free stream velocity ---*/

    if (nDim == 2) {
      config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[1] = sin(Alpha)*Mach*SoundSpeed/Velocity_Ref;
    }
    if (nDim == 3) {
      config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*cos(Beta)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[1] = sin(Beta)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[2] = sin(Alpha)*cos(Beta)*Mach*SoundSpeed/Velocity_Ref;
    }

    config->SetTemperature_FreeStreamND(Temperature);
    direct_solver->SetTemperature_Inf(Temperature);
    config->SetPressure_FreeStreamND(Pressure);
    direct_solver->SetPressure_Inf(Pressure);

  }

  if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS) && config->GetBoolTurbomachinery()){

    BPressure = config->GetPressureOut_BC();
    Temperature = config->GetTotalTemperatureIn_BC();

    if (!reset){
      AD::RegisterInput(BPressure);
      AD::RegisterInput(Temperature);
    }

    config->SetPressureOut_BC(BPressure);
    config->SetTotalTemperatureIn_BC(Temperature);
  }


    /*--- Here it is possible to register other variables as input that influence the flow solution
     * and thereby also the objective function. The adjoint values (i.e. the derivatives) can be
     * extracted in the ExtractAdjointVariables routine. ---*/
}

void CDiscAdjSolver::RegisterOutput(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint, nPoint = geometry->GetnPoint();

  /*--- Register variables as output of the solver iteration ---*/

  bool input = false;

  /*--- Register output variables on the tape ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->node[iPoint]->RegisterSolution(input);
  }
}

void CDiscAdjSolver::RegisterObj_Func(CConfig *config) {

  /*--- Here we can add new (scalar) objective functions ---*/
  if (config->GetnObj()==1) {
    switch (config->GetKind_ObjFunc()) {
    case DRAG_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CD();
      if (config->GetFixed_CL_Mode()) ObjFunc_Value -= config->GetdCD_dCL() * direct_solver->GetTotal_CL();
      if (config->GetFixed_CM_Mode()) ObjFunc_Value -= config->GetdCD_dCMy() * direct_solver->GetTotal_CMy();
      break;
    case LIFT_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CL();
      break;
    case SIDEFORCE_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CSF();
      break;
    case EFFICIENCY:
      ObjFunc_Value = direct_solver->GetTotal_CEff();
      break;
    case MOMENT_X_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMx();
      break;
    case MOMENT_Y_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMy();
      break;
    case MOMENT_Z_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMz();
      break;
    case EQUIVALENT_AREA:
      ObjFunc_Value = direct_solver->GetTotal_CEquivArea();
      break;
    }

    /*--- Template for new objective functions where TemplateObjFunction()
     *  is the routine that returns the obj. function value. The computation
     * must be done while the tape is active, i.e. between AD::StartRecording() and
     * AD::StopRecording() in DiscAdjMeanFlowIteration::Iterate(). The best place is somewhere
     * inside MeanFlowIteration::Iterate().
     *
     * case TEMPLATE_OBJECTIVE:
     *    ObjFunc_Value = TemplateObjFunction();
     *    break;
     * ---*/
  }
  else{
    ObjFunc_Value = direct_solver->GetTotal_ComboObj();
  }
  if (rank == MASTER_NODE) {
    AD::RegisterOutput(ObjFunc_Value);
  }
}

void CDiscAdjSolver::SetAdj_ObjFunc(CGeometry *geometry, CConfig *config) {

  bool time_stepping = config->GetUnsteady_Simulation() != STEADY;
  unsigned long IterAvg_Obj = config->GetIter_Avg_Objective();
  unsigned long ExtIter = config->GetExtIter();
  su2double seeding = 1.0;

  if (time_stepping) {
    if (ExtIter < IterAvg_Obj) {
      seeding = 1.0/((su2double)IterAvg_Obj);
    }
    else {
      seeding = 0.0;
    }
  }

  if (rank == MASTER_NODE) {
    SU2_TYPE::SetDerivative(ObjFunc_Value, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc_Value, 0.0);
  }
}

void CDiscAdjSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){

  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  bool time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
      SetRes_RMS(iVar,0.0);
      SetRes_Max(iVar,0.0,0);
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Set the old solution ---*/

    node[iPoint]->Set_OldSolution();

    /*--- Extract the adjoint solution ---*/

    direct_solver->node[iPoint]->GetAdjointSolution(Solution);

    /*--- Store the adjoint solution ---*/

    node[iPoint]->SetSolution(Solution);
  }

  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {

      /*--- Extract the adjoint solution at time n ---*/

      direct_solver->node[iPoint]->GetAdjointSolution_time_n(Solution);

      /*--- Store the adjoint solution at time n ---*/

      node[iPoint]->Set_Solution_time_n(Solution);
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {

      /*--- Extract the adjoint solution at time n-1 ---*/

      direct_solver->node[iPoint]->GetAdjointSolution_time_n1(Solution);

      /*--- Store the adjoint solution at time n-1 ---*/

      node[iPoint]->Set_Solution_time_n1(Solution);
    }
  }

  /*--- Set the residuals ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
          residual = node[iPoint]->GetSolution(iVar) - node[iPoint]->GetSolution_Old(iVar);

          AddRes_RMS(iVar,residual*residual);
          AddRes_Max(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
  }

  SetResidual_RMS(geometry, config);
}

void CDiscAdjSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) {

  /*--- Extract the adjoint values of the farfield values ---*/

  if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS) && !config->GetBoolTurbomachinery()) {
    su2double Local_Sens_Press, Local_Sens_Temp, Local_Sens_AoA, Local_Sens_Mach;

    Local_Sens_Mach  = SU2_TYPE::GetDerivative(Mach);
    Local_Sens_AoA   = SU2_TYPE::GetDerivative(Alpha);
    Local_Sens_Temp  = SU2_TYPE::GetDerivative(Temperature);
    Local_Sens_Press = SU2_TYPE::GetDerivative(Pressure);

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&Local_Sens_Mach,  &Total_Sens_Mach,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_AoA,   &Total_Sens_AoA,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_Temp,  &Total_Sens_Temp,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    Total_Sens_Mach  = Local_Sens_Mach;
    Total_Sens_AoA   = Local_Sens_AoA;
    Total_Sens_Temp  = Local_Sens_Temp;
    Total_Sens_Press = Local_Sens_Press;
#endif
  }

  if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS) && config->GetBoolTurbomachinery()){
    su2double Local_Sens_BPress, Local_Sens_Temperature;

    Local_Sens_BPress = SU2_TYPE::GetDerivative(BPressure);
    Local_Sens_Temperature = SU2_TYPE::GetDerivative(Temperature);

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&Local_Sens_BPress,   &Total_Sens_BPress,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_Temperature,   &Total_Sens_Temp,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    Total_Sens_BPress = Local_Sens_BPress;
    Total_Sens_Temp = Local_Sens_Temperature;
#endif

  }

  /*--- Extract here the adjoint values of everything else that is registered as input in RegisterInput. ---*/

}


void CDiscAdjSolver::ExtractAdjoint_Geometry(CGeometry *geometry, CConfig *config) {

//  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
//      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

//  bool time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

//  unsigned short iVar;
  unsigned long iPoint;

  /*--- Set Residuals to zero ---*/

//  for (iVar = 0; iVar < nVar; iVar++){
//      SetRes_RMS(iVar,0.0);
//      SetRes_Max(iVar,0.0,0);
//  }

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Set the old solution ---*/

    node[iPoint]->Set_OldSolution_Geometry();

    /*--- Extract the adjoint solution ---*/

    geometry->node[iPoint]->GetAdjointCoord(Solution_Geometry);

    /*--- Store the adjoint solution ---*/

    node[iPoint]->SetSolution_Geometry(Solution_Geometry);

  }

//  if (time_n_needed){
//    for (iPoint = 0; iPoint < nPoint; iPoint++){
//
//      /*--- Extract the adjoint solution at time n ---*/
//
//      direct_solver->node[iPoint]->GetAdjointSolution_time_n(Solution);
//
//      /*--- Store the adjoint solution at time n ---*/
//
//      node[iPoint]->Set_Solution_time_n(Solution);
//    }
//  }
//  if (time_n1_needed){
//    for (iPoint = 0; iPoint < nPoint; iPoint++){
//
//      /*--- Extract the adjoint solution at time n-1 ---*/
//
//      direct_solver->node[iPoint]->GetAdjointSolution_time_n1(Solution);
//
//      /*--- Store the adjoint solution at time n-1 ---*/
//
//      node[iPoint]->Set_Solution_time_n1(Solution);
//    }
//  }

  /*--- Set the residuals ---*/

//  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
//      for (iVar = 0; iVar < nVar; iVar++){
//          residual = node[iPoint]->GetSolution_Geometry(iVar) - node[iPoint]->Get_OldSolution_Geometry(iVar);
//
//          AddRes_RMS(iVar,residual*residual);
//          AddRes_Max(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
//      }
//  }
//
//  SetResidual_RMS(geometry, config);
}

void CDiscAdjSolver::ExtractAdjoint_CrossTerm(CGeometry *geometry, CConfig *config) {

  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution ---*/

    direct_solver->node[iPoint]->GetAdjointSolution(Solution);

    for (iVar = 0; iVar < nVar; iVar++) node[iPoint]->SetCross_Term_Derivative(iVar, Solution[iVar]);

  }

//  Set_MPI_CrossTerm(geometry, config);

}

void CDiscAdjSolver::ExtractAdjoint_CrossTerm_Geometry(CGeometry *geometry, CConfig *config) {

  unsigned short iDim;
  unsigned long iPoint;


  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution ---*/

    geometry->node[iPoint]->GetAdjointCoord(Solution_Geometry);

    for (iDim = 0; iDim < nDim; iDim++) node[iPoint]->SetGeometry_CrossTerm_Derivative(iDim, Solution_Geometry[iDim]);

  }

}

void CDiscAdjSolver::ExtractAdjoint_CrossTerm_Geometry_Flow(CGeometry *geometry, CConfig *config){

  unsigned short iDim;
  unsigned long iPoint;


  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution ---*/

    geometry->node[iPoint]->GetAdjointCoord(Solution_Geometry);

    for (iDim = 0; iDim < nDim; iDim++) node[iPoint]->SetGeometry_CrossTerm_Derivative_Flow(iDim, Solution_Geometry[iDim]);

  }

}


void CDiscAdjSolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config) {

  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST ||
      config->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool fsi = config->GetFSI_Simulation();

  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Solution[iVar] = node[iPoint]->GetSolution(iVar);
    }
    if (fsi) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution[iVar] += node[iPoint]->GetCross_Term_Derivative(iVar);
      }
    }
    if (dual_time) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution[iVar] += node[iPoint]->GetDual_Time_Derivative(iVar);
      }
    }
    direct_solver->node[iPoint]->SetAdjointSolution(Solution);
  }

}

void CDiscAdjSolver::SetAdjoint_OutputMesh(CGeometry *geometry, CConfig *config){

//  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST ||
//      config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  bool fsi = config->GetFSI_Simulation();

  unsigned short iDim;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iDim = 0; iDim < nDim; iDim++){
      Solution_Geometry[iDim] = 0.0;
    }
    if (fsi){
      for (iDim = 0; iDim < nDim; iDim++){
        Solution_Geometry[iDim] += node[iPoint]->GetGeometry_CrossTerm_Derivative(iDim);
      }
      for (iDim = 0; iDim < nDim; iDim++){
        Solution_Geometry[iDim] += node[iPoint]->GetGeometry_CrossTerm_Derivative_Flow(iDim);
      }
    }
//    if (dual_time){
//      for (iDim = 0; iDim < nVar; iDim++){
//        Solution_Geometry[iDim] += node[iPoint]->GetDual_Time_Derivative_Geometry(iDim);
//      }
//    }
    geometry->node[iPoint]->SetAdjointCoord(Solution_Geometry);
  }

}

void CDiscAdjSolver::SetSensitivity(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned short iDim;
  su2double *Coord, Sensitivity, eps;

  bool time_stepping = (config->GetUnsteady_Simulation() != STEADY);

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    Coord = geometry->node[iPoint]->GetCoord();

    for (iDim = 0; iDim < nDim; iDim++) {

      Sensitivity = SU2_TYPE::GetDerivative(Coord[iDim]);

      /*--- Set the index manually to zero. ---*/

     AD::ResetInput(Coord[iDim]);

      /*--- If sharp edge, set the sensitivity to 0 on that region ---*/

      if (config->GetSens_Remove_Sharp()) {
        eps = config->GetVenkat_LimiterCoeff()*config->GetRefElemLength();
        if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetAdjSharp_LimiterCoeff()*eps )
          Sensitivity = 0.0;
      }
      if (!time_stepping) {
        node[iPoint]->SetSensitivity(iDim, Sensitivity);
      } else {
        node[iPoint]->SetSensitivity(iDim, node[iPoint]->GetSensitivity(iDim) + Sensitivity);
      }
    }
  }
  SetSurface_Sensitivity(geometry, config);
}

void CDiscAdjSolver::SetSurface_Sensitivity(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, iDim, iMarker_Monitoring;
  unsigned long iVertex, iPoint;
  su2double *Normal, Prod, Sens = 0.0, SensDim, Area, Sens_Vertex, *Sens_Geo;
  Total_Sens_Geo = 0.0;
  string Monitoring_Tag, Marker_Tag;

  Sens_Geo = new su2double[config->GetnMarker_Monitoring()];
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Sens_Geo[iMarker_Monitoring] = 0.0;
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    /*--- Loop over boundary markers to select those for Euler walls and NS walls ---*/

    if(config->GetMarker_All_KindBC(iMarker) == EULER_WALL
       || config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX
       || config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) {

      Sens = 0.0;

      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Prod = 0.0;
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          /*--- retrieve the gradient calculated with AD -- */
          SensDim = node[iPoint]->GetSensitivity(iDim);

          /*--- calculate scalar product for projection onto the normal vector ---*/
          Prod += Normal[iDim]*SensDim;

          Area += Normal[iDim]*Normal[iDim];
        }

        Area = sqrt(Area);

        /*--- Projection of the gradient calculated with AD onto the normal vector of the surface ---*/

        Sens_Vertex = Prod/Area;
        CSensitivity[iMarker][iVertex] = -Sens_Vertex;
        Sens += Sens_Vertex*Sens_Vertex;
      }

      if (config->GetMarker_All_Monitoring(iMarker) == YES){

        /*--- Compute sensitivity for each surface point ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Sens_Geo[iMarker_Monitoring] = Sens;
          }
        }
      }
    }
  }

#ifdef HAVE_MPI
  su2double *MySens_Geo;
  MySens_Geo = new su2double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySens_Geo[iMarker_Monitoring] = Sens_Geo[iMarker_Monitoring];
    Sens_Geo[iMarker_Monitoring]   = 0.0;
  }

  SU2_MPI::Allreduce(MySens_Geo, Sens_Geo, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  delete [] MySens_Geo;
#endif

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Sens_Geo[iMarker_Monitoring] = sqrt(Sens_Geo[iMarker_Monitoring]);
    Total_Sens_Geo   += Sens_Geo[iMarker_Monitoring];
  }
  
  delete [] Sens_Geo;

}

void CDiscAdjSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config_container, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  bool dual_time_1st = (config_container->GetUnsteady_Simulation() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  su2double *solution_n, *solution_n1;
  unsigned long iPoint;
  unsigned short iVar;
  if (dual_time) {
      for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++) {
          solution_n = node[iPoint]->GetSolution_time_n();
          solution_n1 = node[iPoint]->GetSolution_time_n1();
          for (iVar=0; iVar < nVar; iVar++) {
              node[iPoint]->SetDual_Time_Derivative(iVar, solution_n[iVar]+node[iPoint]->GetDual_Time_Derivative_n(iVar));
              node[iPoint]->SetDual_Time_Derivative_n(iVar, solution_n1[iVar]);
            }
        }
    }
}

void CDiscAdjSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  unsigned short iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine, counter;
  su2double Area_Children, Area_Parent, *Solution_Fine;
  ifstream restart_file;
  string restart_filename, filename, text_line;

  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  /*--- Restart the solution from file information ---*/

  filename = config->GetSolution_AdjFileName();
  restart_filename = config->GetObjFunc_Extension(filename);

  /*--- Read and store the restart metadata. ---*/

  Read_SU2_Restart_Metadata(geometry[MESH_0], config, true, restart_filename);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Read all lines in the restart file ---*/

  long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Skip coordinates ---*/
  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Skip flow adjoint variables ---*/
  if (KindDirect_Solver== RUNTIME_TURB_SYS) {
    if (compressible) {
      skipVars += nDim + 2;
    }
    if (incompressible) {
      skipVars += nDim + 1;
    }
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
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][ADJFLOW_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][ADJFLOW_SOL]->node[iPoint]->SetSolution(Solution);
    }
  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

void CDiscAdjSolver::ComputeResidual_BGS(CGeometry *geometry, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual, bgs_sol;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++){
      SetRes_BGS(iVar,0.0);
      SetRes_Max_BGS(iVar,0.0,0);
  }

  /*--- Compute the BGS solution (adding the cross term) ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      bgs_sol = node[iPoint]->GetSolution(iVar) + node[iPoint]->GetCross_Term_Derivative(iVar);
      node[iPoint]->Set_BGSSolution(iVar, bgs_sol);
    }
  }

  /*--- Set the residuals ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
          residual = node[iPoint]->Get_BGSSolution(iVar) - node[iPoint]->Get_BGSSolution_k(iVar);

          AddRes_BGS(iVar,residual*residual);
          AddRes_Max_BGS(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
  }

  SetResidual_BGS(geometry, config);

}

void CDiscAdjSolver::UpdateSolution_BGS(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  /*--- To nPoint: The solution must be communicated beforehand ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++){

    node[iPoint]->Set_BGSSolution_k();

  }

}
