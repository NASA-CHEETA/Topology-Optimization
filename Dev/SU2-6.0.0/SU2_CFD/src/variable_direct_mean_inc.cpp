/*!
 * \file variable_direct_mean_inc.cpp
 * \brief Definition of the variable classes for incompressible flow.
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

#include "../include/variable_structure.hpp"

CIncEulerVariable::CIncEulerVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
  
  Primitive = NULL;
  
  Gradient_Primitive = NULL;
  
  Limiter_Primitive = NULL;
  
  WindGust    = NULL;
  WindGustDer = NULL;

  nPrimVar     = 0;
  nPrimVarGrad = 0;

  nSecondaryVar     = 0;
  nSecondaryVarGrad = 0;
 
  Undivided_Laplacian = NULL;
 
}

CIncEulerVariable::CIncEulerVariable(su2double val_pressure, su2double *val_velocity, unsigned short val_nDim,
                               unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {
  unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool windgust = config->GetWind_Gust();
  
  /*--- Array initialization ---*/
  
  Primitive = NULL;
  
  Gradient_Primitive = NULL;
  
  Limiter_Primitive = NULL;
  
  WindGust    = NULL;
  WindGustDer = NULL;
  
  nPrimVar     = 0;
  nPrimVarGrad = 0;
  
  nSecondaryVar     = 0;
  nSecondaryVarGrad = 0;

  Undivided_Laplacian = NULL;

  /*--- Allocate and initialize the primitive variables and gradients ---*/
  
  nPrimVar = nDim+5; nPrimVarGrad = nDim+3;

  /*--- Allocate residual structures ---*/
  
  Res_TruncError = new su2double [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }
  
  /*--- Only for residual smoothing (multigrid) ---*/
  
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }
  
  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    Undivided_Laplacian = new su2double [nVar];
  }
  
  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/

  Limiter_Primitive = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Limiter_Primitive[iVar] = 0.0;

  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;
  
  Solution_Max = new su2double [nPrimVarGrad];
  Solution_Min = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }
  
  /*--- Solution and old solution initialization ---*/

  Solution[0] = val_pressure;
  Solution_Old[0] = val_pressure;
  for (iDim = 0; iDim < nDim; iDim++) {
    Solution[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
    Solution_Old[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
  }

  
  /*--- Allocate and initialize solution for dual time strategy ---*/
  
  if (dual_time) {
    Solution_time_n[0]  =  val_pressure;
    Solution_time_n1[0] =  val_pressure;
    for (iDim = 0; iDim < nDim; iDim++) {
      Solution_time_n[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
      Solution_time_n1[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
    }
  }
    
  /*--- Allocate vector for wind gust and wind gust derivative field ---*/
  
  if (windgust) {
    WindGust = new su2double [nDim];
    WindGustDer = new su2double [nDim+1];
  }

  /*--- Incompressible flow, primitive variables nDim+3, (P, vx, vy, vz, rho, beta) ---*/
  
  Primitive = new su2double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;

  /*--- Incompressible flow, gradients primitive variables nDim+2, (P, vx, vy, vz, rho)
        We need P, and rho for running the adjoint problem ---*/
  
  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }
}

CIncEulerVariable::CIncEulerVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {
  unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool windgust = config->GetWind_Gust();
  
  /*--- Array initialization ---*/
  
  Primitive = NULL;
  
  Gradient_Primitive = NULL;
  
  Limiter_Primitive = NULL;
  
  WindGust    = NULL;
  WindGustDer = NULL;
  
  nPrimVar     = 0;
  nPrimVarGrad = 0;
  
  nSecondaryVar     = 0;
  nSecondaryVarGrad = 0;
 
  Undivided_Laplacian = NULL;
 
  /*--- Allocate and initialize the primitive variables and gradients ---*/
  nPrimVar = nDim+5; nPrimVarGrad = nDim+3;
  
  /*--- Allocate residual structures ---*/
  Res_TruncError = new su2double [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }
  
  /*--- Only for residual smoothing (multigrid) ---*/
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }
  
  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    Undivided_Laplacian = new su2double [nVar];
  
  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/
  
  Limiter_Primitive = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Limiter_Primitive[iVar] = 0.0;

  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;
  
  Solution_Max = new su2double [nPrimVarGrad];
  Solution_Min = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }
  
  /*--- Solution initialization ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar] = val_solution[iVar];
    Solution_Old[iVar] = val_solution[iVar];
  }
  
  /*--- Allocate and initializate solution for dual time strategy ---*/
  
  if (dual_time) {
    Solution_time_n = new su2double [nVar];
    Solution_time_n1 = new su2double [nVar];
    
    for (iVar = 0; iVar < nVar; iVar++) {
      Solution_time_n[iVar] = val_solution[iVar];
      Solution_time_n1[iVar] = val_solution[iVar];
    }
  }
  
    
  /*--- Allocate vector for wind gust and wind gust derivative field ---*/
  
  if (windgust) {
    WindGust = new su2double [nDim];
    WindGustDer = new su2double [nDim+1];
  }
  
  /*--- Incompressible flow, primitive variables nDim+3, (P, vx, vy, vz, rho, beta) ---*/
  
  Primitive = new su2double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;

  /*--- Incompressible flow, gradients primitive variables nDim+2, (P, vx, vy, vz, rho),
        We need P, and rho for running the adjoint problem ---*/
  
  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }
  
}

CIncEulerVariable::~CIncEulerVariable(void) {
  unsigned short iVar;

  if (Primitive         != NULL) delete [] Primitive;
  if (Limiter_Primitive != NULL) delete [] Limiter_Primitive;
  if (WindGust          != NULL) delete [] WindGust;
  if (WindGustDer       != NULL) delete [] WindGustDer;

  if (Gradient_Primitive != NULL) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      if (Gradient_Primitive!=NULL) delete [] Gradient_Primitive[iVar];
    delete [] Gradient_Primitive;
  }

  if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
  
}

void CIncEulerVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) {
  unsigned short iVar, iDim;
  
  for (iVar = 0; iVar < val_primvar; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
}


su2double CIncEulerVariable::GetProjVel(su2double *val_vector) {
  su2double ProjVel;
  unsigned short iDim;
  
  ProjVel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ProjVel += Primitive[iDim+1]*val_vector[iDim];
  
  return ProjVel;
}


bool CIncEulerVariable::SetPrimVar(su2double Density_Inf, CConfig *config) {
  
  su2double ArtComp_Factor = config->GetArtComp_Factor();
  
  /*--- Set the value of the density ---*/
  
  SetDensity(Density_Inf);
  
  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/
  
  SetVelocity();
  
  /*--- Set the value of the pressure ---*/
  
  SetPressure();
  
  /*--- Set the value of the artificial compressibility factor ---*/
  
  SetBetaInc2(ArtComp_Factor);
  
  return true;
  
}

CIncNSVariable::CIncNSVariable(void) : CIncEulerVariable() { }

CIncNSVariable::CIncNSVariable(su2double val_pressure, su2double *val_velocity,
                         unsigned short val_nDim, unsigned short val_nvar,
                         CConfig *config) : CIncEulerVariable(val_pressure, val_velocity, val_nDim, val_nvar, config) {
  
  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();
  
  DES_LengthScale = 0.0;
  
}

CIncNSVariable::CIncNSVariable(su2double *val_solution, unsigned short val_nDim,
                         unsigned short val_nvar, CConfig *config) : CIncEulerVariable(val_solution, val_nDim, val_nvar, config) {
  
  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();
  
  DES_LengthScale = 0.0;
  
}

CIncNSVariable::~CIncNSVariable(void) { }

bool CIncNSVariable::SetVorticity(void) {
  
  Vorticity[0] = 0.0; Vorticity[1] = 0.0;
  
  Vorticity[2] = Gradient_Primitive[2][0]-Gradient_Primitive[1][1];
  
  if (nDim == 3) {
    Vorticity[0] = Gradient_Primitive[3][1]-Gradient_Primitive[2][2];
    Vorticity[1] = -(Gradient_Primitive[3][0]-Gradient_Primitive[1][2]);
  }
  
  return false;
  
}

bool CIncNSVariable::SetStrainMag(void) {
  
  su2double Div;
  unsigned short iDim;
  
  AD::StartPreacc();
  AD::SetPreaccIn(Gradient_Primitive, nDim+1, nDim);

  Div = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Div += Gradient_Primitive[iDim+1][iDim];
  }
  
  StrainMag = 0.0;
  
  /*--- Add diagonal part ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    StrainMag += pow(Gradient_Primitive[iDim+1][iDim] - 1.0/3.0*Div, 2.0);
  }
  
  /*--- Add off diagonals ---*/

  StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][1] + Gradient_Primitive[2][0]), 2.0);

  if (nDim == 3) {
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][2] + Gradient_Primitive[3][0]), 2.0);
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[2][2] + Gradient_Primitive[3][1]), 2.0);
  }
  
  StrainMag = sqrt(2.0*StrainMag);

  AD::SetPreaccOut(StrainMag);
  AD::EndPreacc();

  return false;
  
}


bool CIncNSVariable::SetPrimVar(su2double Density_Inf, su2double Viscosity_Inf, su2double eddy_visc, su2double turb_ke, CConfig *config) {
  
  su2double ArtComp_Factor = config->GetArtComp_Factor();
  
  /*--- Set the value of the density and viscosity ---*/
  
  SetDensity(Density_Inf);
  SetLaminarViscosity(Viscosity_Inf);
  
  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/
  
  SetVelocity();
  
  /*--- Set the value of the pressure ---*/
  
  SetPressure();
  
  /*--- Set the value of the artificial compressibility factor ---*/
  
  SetBetaInc2(ArtComp_Factor);
  
  /*--- Set eddy viscosity ---*/
  
  SetEddyViscosity(eddy_visc);
  
  return true;
  
}

