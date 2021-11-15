/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems
 * \author F. Palacios, A. Bueno
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

CTurbSolver::CTurbSolver(void) : CSolver() {
  
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  lowerlimit    = NULL;
  upperlimit    = NULL;
  nVertex       = NULL;
  nMarker       = 0;
  Inlet_TurbVars = NULL;
  
}

CTurbSolver::CTurbSolver(CGeometry* geometry, CConfig *config) : CSolver() {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  lowerlimit    = NULL;
  upperlimit    = NULL;
  nMarker       = config->GetnMarker_All();
  
  /*--- Store the number of vertices on each marker for deallocation later ---*/
  nVertex = new unsigned long[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];
  
}

CTurbSolver::~CTurbSolver(void) {

  if (Inlet_TurbVars != NULL) {
    for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Inlet_TurbVars[iMarker] != NULL) {
        for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          delete [] Inlet_TurbVars[iMarker][iVertex];
        }
        delete [] Inlet_TurbVars[iMarker];
      }
    }
    delete [] Inlet_TurbVars;
  }

  if (FlowPrimVar_i != NULL) delete [] FlowPrimVar_i;
  if (FlowPrimVar_j != NULL) delete [] FlowPrimVar_j;
  if (lowerlimit != NULL) delete [] lowerlimit;
  if (upperlimit != NULL) delete [] upperlimit;
  if (nVertex != NULL) delete [] nVertex;
  

}

void CTurbSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, nBufferS_Scalar, nBufferR_Scalar;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL, *Buffer_Receive_muT = NULL, *Buffer_Send_muT = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
  
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      nBufferS_Scalar = nVertexS;             nBufferR_Scalar = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      Buffer_Receive_muT = new su2double [nBufferR_Scalar];
      Buffer_Send_muT = new su2double[nBufferS_Scalar];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_muT[iVertex] = node[iPoint]->GetmuT();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      SU2_MPI::Sendrecv(Buffer_Send_muT, nBufferS_Scalar, MPI_DOUBLE, send_to, 1,
                   Buffer_Receive_muT, nBufferR_Scalar, MPI_DOUBLE, receive_from, 1, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        Buffer_Receive_muT[iVertex] = node[iPoint]->GetmuT();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      delete [] Buffer_Send_muT;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy conservative variables. ---*/
        node[iPoint]->SetmuT(Buffer_Receive_muT[iVertex]);
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_muT;
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CTurbSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
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
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
}

void CTurbSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
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
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
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

void CTurbSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nVar];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
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
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetLimiter(iVar, Buffer_Receive_Limit[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
  }
  
  delete [] Limiter;
  
}


void CTurbSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
  
  su2double *Turb_i, *Turb_j, *Limiter_i = NULL, *Limiter_j = NULL, *V_i, *V_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;
  
  bool muscl         = config->GetMUSCL_Turb();
  bool limiter       = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER);
  bool grid_movement = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive variables w/o reconstruction ---*/
    
    V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
    V_j = solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive();
    numerics->SetPrimitive(V_i, V_j);
    
    /*--- Turbulent variables w/o reconstruction ---*/
    
    Turb_i = node[iPoint]->GetSolution();
    Turb_j = node[jPoint]->GetSolution();
    numerics->SetTurbVar(Turb_i, Turb_j);
    
    /*--- Grid Movement ---*/
    
    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    if (muscl) {

      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      /*--- Mean flow primitive variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
      Gradient_j = solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive();
      if (limiter) {
        Limiter_i = solver_container[FLOW_SOL]->node[iPoint]->GetLimiter_Primitive();
        Limiter_j = solver_container[FLOW_SOL]->node[jPoint]->GetLimiter_Primitive();
      }
      
      for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnPrimVarGrad(); iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          FlowPrimVar_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          FlowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetPrimitive(FlowPrimVar_i, FlowPrimVar_j);
      
      /*--- Turbulent variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = node[iPoint]->GetGradient();
      Gradient_j = node[jPoint]->GetGradient();
      if (limiter) {
        Limiter_i = node[iPoint]->GetLimiter();
        Limiter_j = node[jPoint]->GetLimiter();
      }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          Solution_i[iVar] = Turb_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          Solution_j[iVar] = Turb_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Solution_i[iVar] = Turb_i[iVar] + Project_Grad_i;
          Solution_j[iVar] = Turb_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetTurbVar(Solution_i, Solution_j);
      
    }
    
    /*--- Add and subtract residual ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.AddBlock(iPoint, Residual);
    LinSysRes.SubtractBlock(jPoint, Residual);
    
    /*--- Implicit part ---*/
    
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CTurbSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iEdge, iPoint, jPoint;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and normal vector ---*/
    
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(),
                           solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive());
    
    /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
    
    numerics->SetTurbVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
    /*--- Menter's first blending function (only SST)---*/
    if (config->GetKind_Turb_Model() == SST)
      numerics->SetF1blending(node[iPoint]->GetF1blending(), node[jPoint]->GetF1blending());
    
    /*--- Compute residual, and Jacobians ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    /*--- Add and subtract residual, and update Jacobians ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CTurbSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  /*--- Convective fluxes across symmetry plane are equal to zero. ---*/

}

void CTurbSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  
  /*--- Convective fluxes across euler wall are equal to zero. ---*/

}

void CTurbSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Riemann(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
    BC_Inlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case STATIC_PRESSURE:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}

void CTurbSolver::BC_TurboRiemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Riemann(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
    BC_Inlet_Turbo(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case STATIC_PRESSURE:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}


void CTurbSolver::BC_Giles(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Giles(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT:case TOTAL_CONDITIONS_PT_1D: case DENSITY_VELOCITY:
    BC_Inlet_Turbo(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case MIXING_IN:
    if (config->GetBoolTurbMixingPlane()){
      BC_Inlet_MixingPlane(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    }
    else{
      BC_Inlet_Turbo(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    }
    break;

  case STATIC_PRESSURE: case MIXING_OUT: case STATIC_PRESSURE_1D: case RADIAL_EQUILIBRIUM:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}

void CTurbSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, Vol, density_old = 0.0, density = 0.0;
  
  bool adjoint = config->GetContinuous_Adjoint() || (config->GetDiscrete_Adjoint() && config->GetFrozen_Visc_Disc());
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the volume ---*/
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    Delta = Vol / (config->GetCFLRedCoeff_Turb()*solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    Jacobian.AddVal2Diag(iPoint, Delta);
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - LinSysRes[total_index];
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
  system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  if (!adjoint) {
    
    /*--- Update and clip trubulent solution ---*/
    
    switch (config->GetKind_Turb_Model()) {
        
      case SA: case SA_E: case SA_COMP: case SA_E_COMP: 
        
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          node[iPoint]->AddClippedSolution(0, config->GetRelaxation_Factor_Turb()*LinSysSol[iPoint], lowerlimit[0], upperlimit[0]);
        }
        
        break;
        
      case SA_NEG:
        
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          node[iPoint]->AddSolution(0, config->GetRelaxation_Factor_Turb()*LinSysSol[iPoint]);
        }
        
        break;

      case SST:
        
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          
          if (compressible) {
            density_old = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_Old(0);
            density     = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
          }
          if (incompressible) {
            density_old = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
            density     = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
          }
          
          for (iVar = 0; iVar < nVar; iVar++) {
            node[iPoint]->AddConservativeSolution(iVar, config->GetRelaxation_Factor_Turb()*LinSysSol[iPoint*nVar+iVar], density, density_old, lowerlimit[iVar], upperlimit[iVar]);
          }
          
        }
        
        break;
        
    }
  }
  
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CTurbSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                       unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;
  
  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
  su2double Volume_nM1, Volume_nP1, TimeStep;
  su2double Density_nM1, Density_n, Density_nP1;
  su2double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  
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
      
      if (config->GetKind_Turb_Model() == SST) {
        
        /*--- If this is the SST model, we need to multiply by the density
         in order to get the conservative variables ---*/
        Density_nM1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[0];
        Density_n   = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
        Density_nP1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution()[0];
        
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Residual[iVar] = ( Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Residual[iVar] = ( 3.0*Density_nP1*U_time_nP1[iVar] - 4.0*Density_n*U_time_n[iVar]
                              +1.0*Density_nM1*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
        }
        
      } else {
        
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                              +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
        }
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
    
  } else {
    
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
      
      /*--- Multiply by density at node i for the SST model ---*/
      
      if (config->GetKind_Turb_Model() == SST) {
        Density_n = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;
      } else {
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      }
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Compute the GCL component of the source term for node j ---*/
      
      U_time_n = node[jPoint]->GetSolution_time_n();
      
      /*--- Multiply by density at node j for the SST model ---*/
      
      if (config->GetKind_Turb_Model() == SST) {
        Density_n = solver_container[FLOW_SOL]->node[jPoint]->GetSolution_time_n()[0];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;
      } else {
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      }
      LinSysRes.SubtractBlock(jPoint, Residual);
      
    }
    
    /*---  Loop over the boundary edges ---*/
    
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
        
        /*--- Multiply by density at node i for the SST model ---*/
        
        if (config->GetKind_Turb_Model() == SST) {
          Density_n = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;
        } else {
          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = U_time_n[iVar]*Residual_GCL;
        }
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
      
      if (config->GetKind_Turb_Model() == SST) {
        
        /*--- If this is the SST model, we need to multiply by the density
         in order to get the conservative variables ---*/
        Density_nM1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[0];
        Density_n   = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
        Density_nP1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution()[0];
        
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Residual[iVar] = (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(Volume_nP1/TimeStep);
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Residual[iVar] = (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
            + (Density_nM1*U_time_nM1[iVar] - Density_n*U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
        }
        
      } else {
        
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
            + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
        }
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


void CTurbSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  
  unsigned short iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  su2double Area_Children, Area_Parent, *Solution_Fine;
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  string UnstExt, text_line;
  ifstream restart_file;
  string restart_filename = config->GetSolution_FlowFileName();

  /*--- Modify file name for multizone problems ---*/
  if (nZone >1)
    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/
  
  if (dual_time|| time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Skip flow variables ---*/
  
  unsigned short skipVars = 0;

  if (compressible) {
    if (nDim == 2) skipVars += 6;
    if (nDim == 3) skipVars += 8;
  }
  if (incompressible) {
    if (nDim == 2) skipVars += 5;
    if (nDim == 3) skipVars += 7;
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
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- MPI solution and compute the eddy viscosity ---*/

//TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart.
  solver[MESH_0][TURB_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][TURB_SOL]->Set_MPI_Solution(geometry[MESH_0], config);

  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/
  
  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][TURB_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][TURB_SOL]->node[iPoint]->SetSolution(Solution);
    }
    solver[iMesh][TURB_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    solver[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

CTurbSASolver::CTurbSASolver(void) : CTurbSolver() {

  Inlet_TurbVars = NULL;

}

CTurbSASolver::CTurbSASolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, CFluidModel* FluidModel)
    : CTurbSolver(geometry, config) {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  su2double Density_Inf, Viscosity_Inf, Factor_nu_Inf, Factor_nu_Engine, Factor_nu_ActDisk;

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Dimension of the problem --> dependent of the turbulent model ---*/
  
  nVar = 1;
  nPrimVar = 1;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  
  /*--- Single grid simulation ---*/
  
  if (iMesh == MESH_0 || config->GetMGCycle() == FULLMG_CYCLE) {
    
    /*--- Define some auxiliar vector related with the residual ---*/
    
    Residual = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]  = 0.0;
    Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_i = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]  = 0.0;
    Residual_j = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]  = 0.0;
    Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
    
    /*--- Define some structures for locating max residuals ---*/
    
    Point_Max = new unsigned long[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
    Point_Max_Coord = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
    }
    
    /*--- Define some auxiliar vector related with the solution ---*/
    
    Solution = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
    
    /*--- Define some auxiliar vector related with the geometry ---*/
    
    Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
    
    /*--- Define some auxiliar vector related with the flow solution ---*/
    
    FlowPrimVar_i = new su2double [nDim+7]; FlowPrimVar_j = new su2double [nDim+7];
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    /*--- Initialization of the structure of the whole Jacobian ---*/
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SA model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    
    if (config->GetExtraOutput()) {
      if (nDim == 2) { nOutputVariables = 13; }
      else if (nDim == 3) { nOutputVariables = 19; }
      OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
      OutputHeadingNames = new string[nOutputVariables];
    }
    
    /*--- Computation of gradients by least squares ---*/
    
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
      Smatrix = new su2double* [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Smatrix[iDim] = new su2double [nDim];
      
      /*--- c vector := transpose(WA)*(Wb) ---*/
      
      Cvector = new su2double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        Cvector[iVar] = new su2double [nDim];
    }
    
  }
  
  /*--- Initialize lower and upper limits---*/
  
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  
  lowerlimit[0] = 1.0e-10;
  upperlimit[0] = 1.0;
  

  /*--- Read farfield conditions from config ---*/
  
  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  
  /*--- Factor_nu_Inf in [3.0, 5.0] ---*/

  Factor_nu_Inf = config->GetNuFactor_FreeStream();
  nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  if (config->GetKind_Trans_Model() == BC) {
    nu_tilde_Inf  = 0.005*Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  }

  /*--- Factor_nu_Engine ---*/
  Factor_nu_Engine = config->GetNuFactor_Engine();
  nu_tilde_Engine  = Factor_nu_Engine*Viscosity_Inf/Density_Inf;
  if (config->GetKind_Trans_Model() == BC) {
    nu_tilde_Engine  = 0.005*Factor_nu_Engine*Viscosity_Inf/Density_Inf;
  }

  /*--- Factor_nu_ActDisk ---*/
  Factor_nu_ActDisk = config->GetNuFactor_Engine();
  nu_tilde_ActDisk  = Factor_nu_ActDisk*Viscosity_Inf/Density_Inf;

  /*--- Eddy viscosity at infinity ---*/
  su2double Ji, Ji_3, fv1, cv1_3 = 7.1*7.1*7.1;
  su2double muT_Inf;
  Ji = nu_tilde_Inf/Viscosity_Inf*Density_Inf;
  Ji_3 = Ji*Ji*Ji;
  fv1 = Ji_3/(Ji_3+cv1_3);
  muT_Inf = Density_Inf*fv1*nu_tilde_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CTurbSAVariable(nu_tilde_Inf, muT_Inf, nDim, nVar, config);

  /*--- MPI solution ---*/

//TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart
  Set_MPI_Solution(geometry, config);
  Set_MPI_Solution(geometry, config);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  unsigned long iMarker;

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

  /*--- Set up inlet profiles, if necessary ---*/

  SetInlet(config);

}

CTurbSASolver::~CTurbSASolver(void) {
  
  unsigned long iMarker, iVertex;
  unsigned short iVar;
  
  if ( SlidingState != NULL ) {
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

}

void CTurbSASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint;
  unsigned long ExtIter = config->GetExtIter();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter_flow     = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_turb     = ((config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  unsigned short kind_hybridRANSLES = config->GetKind_HybridRANSLES();
  su2double** PrimGrad_Flow = NULL;
  su2double* Vorticity = NULL;
  su2double Laminar_Viscosity = 0;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the residual vector ---*/
    
    LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  Jacobian.SetValZero();

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

  /*--- Upwind second order reconstruction ---*/

  if (limiter_turb) SetSolution_Limiter(geometry, config);

  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);

  if (kind_hybridRANSLES != NO_HYBRIDRANSLES){
    
    /*--- Set the vortex tilting coefficient at every node if required ---*/
    
    if (kind_hybridRANSLES == SA_EDDES){
      for (iPoint = 0; iPoint < nPoint; iPoint++){
        PrimGrad_Flow      = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
        Vorticity          = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity();
        Laminar_Viscosity  = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
        node[iPoint]->SetVortex_Tilting(PrimGrad_Flow, Vorticity, Laminar_Viscosity);
      }
    }
    
    /*--- Compute the DES length scale ---*/
    
    SetDES_LengthScale(solver_container, geometry, config);
    
  }
}

void CTurbSASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  
  su2double rho = 0.0, mu = 0.0, nu, *nu_hat, muT, Ji, Ji_3, fv1;
  su2double cv1_3 = 7.1*7.1*7.1;
  unsigned long iPoint;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool neg_spalart_allmaras = (config->GetKind_Turb_Model() == SA_NEG);
  
  
  /*--- Compute eddy viscosity ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    if (compressible) {
      rho = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      mu  = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    }
    if (incompressible) {
      rho = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      mu  = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    }
    
    nu  = mu/rho;
    nu_hat = node[iPoint]->GetSolution();
    
    Ji   = nu_hat[0]/nu;
    Ji_3 = Ji*Ji*Ji;
    fv1  = Ji_3/(Ji_3+cv1_3);
    
    muT = rho*fv1*nu_hat[0];
    
    if (neg_spalart_allmaras && (muT < 0.0)) muT = 0.0;
    
    node[iPoint]->SetmuT(muT);
    
  }
  
}

void CTurbSASolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                    CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
  
  bool harmonic_balance = (config->GetUnsteady_Simulation() == HARMONIC_BALANCE);
  bool transition    = (config->GetKind_Trans_Model() == LM);
  bool transition_BC = (config->GetKind_Trans_Model() == BC);
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);
    
    /*--- Gradient of the primitive and conservative variables ---*/
    
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
    
    /*--- Set vorticity and strain rate magnitude ---*/
    
    numerics->SetVorticity(solver_container[FLOW_SOL]->node[iPoint]->GetVorticity(), NULL);

    numerics->SetStrainMag(solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag(), 0.0);
    
    /*--- Set intermittency ---*/
    
    if (transition) {
      numerics->SetIntermittency(solver_container[TRANS_SOL]->node[iPoint]->GetIntermittency());
    }
    
    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    
    numerics->SetTurbVar(node[iPoint]->GetSolution(), NULL);
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);
    
    /*--- Set volume ---*/
    
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    /*--- Get Hybrid RANS/LES Type and set the appropriate wall distance ---*/     
    
    if (config->GetKind_HybridRANSLES() == NO_HYBRIDRANSLES) {
          
      /*--- Set distance to the surface ---*/
          
      numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
    
    } else {
    
      /*--- Set DES length scale ---*/
      
      numerics->SetDistance(node[iPoint]->GetDES_LengthScale(), 0.0);
      
    }

    /*--- Compute the source term ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);

    /*--- Store the intermittency ---*/

    if (transition_BC) {
      node[iPoint]->SetGammaBC(numerics->GetGammaBC());
    }
    
    /*--- Subtract residual and the Jacobian ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    
  }
  
  if (harmonic_balance) {
    
    su2double Volume, Source;
    unsigned short nVar_Turb = solver_container[TURB_SOL]->GetnVar();
    
    /*--- Loop over points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Get control volume ---*/
      
      Volume = geometry->node[iPoint]->GetVolume();
      
      /*--- Access stored harmonic balance source term ---*/
      
      for (unsigned short iVar = 0; iVar < nVar_Turb; iVar++) {
        Source = node[iPoint]->GetHarmonicBalance_Source(iVar);
        Residual[iVar] = Source*Volume;
      }
      
      /*--- Add Residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
  }
  
}

void CTurbSASolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                    CConfig *config, unsigned short iMesh) {
  
}

void CTurbSASolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Get the velocity vector ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] = 0.0;
      
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- includes 1 in the diagonal ---*/
      
      Jacobian.DeleteValsRowi(iPoint);
    }
  }
  
}

void CTurbSASolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                       unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Get the velocity vector ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] = 0.0;
      
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Includes 1 in the diagonal ---*/
      
      Jacobian.DeleteValsRowi(iPoint);
    }
  }
  
}

void CTurbSASolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *Normal, *V_infty, *V_domain;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set turbulent variable at the wall, and at infinity ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      Solution_j[0] = nu_tilde_Inf;
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute residuals and Jacobians ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  delete [] Normal;
  
}

void CTurbSASolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_inlet, *V_domain, *Normal;
  bool isCustomizable = config->GetMarker_All_PyCustom(val_marker);
  
  Normal = new su2double[nDim];
  
  bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the inlet ---*/
      
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/
      
      Solution_i[0] = node[iPoint]->GetSolution(0);

      if (isCustomizable) {
        /*--- Only use the array version of the inlet BC when necessary ---*/
        Solution_j[0] = Inlet_TurbVars[val_marker][iVertex][0];
      } else {
        Solution_j[0] = nu_tilde_Inf;
      }
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the conv_numerics class ---*/
      
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_inlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

void CTurbSASolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                              CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the outlet ---*/
      
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TurbVar_internal,
       Solution_j --> TurbVar_outlet ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete[] Normal;
  
}

void CTurbSASolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iDim;
  su2double *V_inflow, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();

  Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      
      V_inflow = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inflow);
      
      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual. ---*/
      
      conv_numerics->SetTurbVar(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
      /*--- Set Normal (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Set grid movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_inflow);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetTurbVar(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }

  }
  
  /*--- Free locally allocated memory ---*/
  
  delete[] Normal;
  
}

void CTurbSASolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_exhaust, *V_domain, *Normal;
  
  Normal = new su2double[nDim];
  
  bool grid_movement  = config->GetGrid_Movement();

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the infinity ---*/
      
      V_exhaust = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_exhaust);
      
      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/
      
      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_j[0] = nu_tilde_Engine;
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the conv_numerics class ---*/
      
      conv_numerics->SetNormal(Normal);

      /*--- Set grid movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_exhaust);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete[] Normal;
  
}

void CTurbSASolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics,
             config,  val_marker, true);
  
}

void CTurbSASolver::BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics,
             config,  val_marker, false);
  
}

void CTurbSASolver::BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                               CConfig *config, unsigned short val_marker, bool inlet_surface) {
  
  unsigned long iPoint, iVertex, GlobalIndex_donor, GlobalIndex;
  su2double *V_outlet, *V_inlet, *V_domain, *Normal, *UnitNormal, Area, Vn;
  bool ReverseFlow;
  unsigned short iDim;
  
  bool grid_movement = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  UnitNormal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex_donor = solver_container[FLOW_SOL]->GetDonorGlobalIndex(val_marker, iVertex);
    GlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex != GlobalIndex_donor)) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Check the flow direction. Project the flow into the normal to the inlet face ---*/
      
      Vn = 0.0; ReverseFlow = false;
      for (iDim = 0; iDim < nDim; iDim++) {  Vn += V_domain[iDim+1]*UnitNormal[iDim]; }
      
      if ((inlet_surface) && (Vn < 0.0)) { ReverseFlow = true; }
      if ((!inlet_surface) && (Vn > 0.0)) { ReverseFlow = true; }
      
      /*--- Do not anything if there is a
       reverse flow, Euler b.c. for the direct problem ---*/
      
      if (!ReverseFlow) {
        
        /*--- Allocate the value at the infinity ---*/
        
        if (inlet_surface) {
          V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
          V_outlet = solver_container[FLOW_SOL]->GetDonorPrimVar(val_marker, iVertex);
          conv_numerics->SetPrimitive(V_domain, V_inlet);
        }
        else {
          V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
          V_inlet = solver_container[FLOW_SOL]->GetDonorPrimVar(val_marker, iVertex);
          conv_numerics->SetPrimitive(V_domain, V_outlet);
        }
        
        /*--- Set the turb. variable solution
         set  the turbulent variables. Here we use a Neumann BC such
         that the turbulent variable is copied from the interior of the
         domain to the outlet before computing the residual.
         or set the turbulent variable states (prescribed for an inflow)  ----*/
        
        Solution_i[0] = node[iPoint]->GetSolution(0);
        
        //      if (inlet_surface) Solution_j[0] = 0.5*(node[iPoint]->GetSolution(0)+V_outlet [nDim+9]);
        //      else Solution_j[0] = 0.5*(node[iPoint]->GetSolution(0)+V_inlet [nDim+9]);
        
        //      /*--- Inflow analysis (interior extrapolation) ---*/
        //      if (((inlet_surface) && (!ReverseFlow)) || ((!inlet_surface) && (ReverseFlow))) {
        //        Solution_j[0] = 2.0*node[iPoint]->GetSolution(0) - node[iPoint_Normal]->GetSolution(0);
        //      }
        
        //      /*--- Outflow analysis ---*/
        //      else {
        //        if (inlet_surface) Solution_j[0] = Factor_nu_ActDisk*V_outlet [nDim+9];
        //        else { Solution_j[0] = Factor_nu_ActDisk*V_inlet [nDim+9]; }
        //      }
        
        /*--- Inflow analysis (interior extrapolation) ---*/
        if (((inlet_surface) && (!ReverseFlow)) || ((!inlet_surface) && (ReverseFlow))) {
          Solution_j[0] = node[iPoint]->GetSolution(0);
        }
        
        /*--- Outflow analysis ---*/
        else {
          Solution_j[0] = nu_tilde_ActDisk;
        }
        
        conv_numerics->SetTurbVar(Solution_i, Solution_j);
        
        /*--- Grid Movement ---*/
        
        if (grid_movement)
          conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
        
        /*--- Compute the residual using an upwind scheme ---*/
        
        conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.AddBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        
//        /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//        visc_numerics->SetNormal(Normal);
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint_Normal]->GetCoord());
//
//        /*--- Conservative variables w/o reconstruction ---*/
//
//        if (inlet_surface) visc_numerics->SetPrimitive(V_domain, V_inlet);
//        else visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//        /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//        visc_numerics->SetTurbVar(Solution_i, Solution_j);
//
//        visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//        /*--- Compute residual, and Jacobians ---*/
//
//        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//        /*--- Subtract residual, and update Jacobians ---*/
//
//        LinSysRes.SubtractBlock(iPoint, Residual);
//        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete[] Normal;
  delete[] UnitNormal;
  
}

void CTurbSASolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim, iSpan;
  unsigned long  oldVertex, iPoint, Point_Normal;
  long iVertex;
  su2double *V_inlet, *V_domain, *Normal;
  su2double extAverageNu;
  Normal = new su2double[nDim];

  bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  unsigned short nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){
    extAverageNu = solver_container[FLOW_SOL]->GetExtAverageNu(val_marker, iSpan);

    /*--- Loop over all the vertices on this boundary marker ---*/

    for (iVertex = 0; iVertex < geometry->nVertexSpan[val_marker][iSpan]; iVertex++) {

      /*--- find the node related to the vertex ---*/
      iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_j[0] = extAverageNu;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
            geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/

      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/

      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/

      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

      /*--- Compute residual, and Jacobians ---*/

      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Subtract residual, and update Jacobians ---*/

      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}

void CTurbSASolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim, iSpan;
  unsigned long  oldVertex, iPoint, Point_Normal;
  long iVertex;
  su2double *V_inlet, *V_domain, *Normal;

  su2double rho, pressure, muLam, Factor_nu_Inf, nu_tilde;
  Normal = new su2double[nDim];

  bool grid_movement  = config->GetGrid_Movement();
  unsigned short nSpanWiseSections = config->GetnSpanWiseSections();
  CFluidModel *FluidModel;

  FluidModel = solver_container[FLOW_SOL]->GetFluidModel();
  Factor_nu_Inf = config->GetNuFactor_FreeStream();


  /*--- Loop over all the spans on this boundary marker ---*/
  for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){
    rho       = solver_container[FLOW_SOL]->GetAverageDensity(val_marker, iSpan);
    pressure  = solver_container[FLOW_SOL]->GetAveragePressure(val_marker, iSpan);

    FluidModel->SetTDState_Prho(pressure, rho);
    muLam = FluidModel->GetLaminarViscosity();

    nu_tilde  = Factor_nu_Inf*muLam/rho;


    /*--- Loop over all the vertices on this boundary marker ---*/
    for (iVertex = 0; iVertex < geometry->nVertexSpan[val_marker][iSpan]; iVertex++) {

      /*--- find the node related to the vertex ---*/
      iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_j[0] =  nu_tilde;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
            geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/

      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/

      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/

      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

      /*--- Compute residual, and Jacobians ---*/

      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Subtract residual, and update Jacobians ---*/

      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}

void CTurbSASolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                          CConfig *config, unsigned short val_marker) {
  
  //  unsigned long iVertex, iPoint, jPoint;
  //  unsigned short iVar, iDim;
  //
  //  su2double *Vector = new su2double[nDim];
  //
  //#ifndef HAVE_MPI
  //
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();
  //
  //      if (iPoint != jPoint) {
  //
  //        /*--- Store the solution for both points ---*/
  //        for (iVar = 0; iVar < nVar; iVar++) {
  //          Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
  //          Solution_j[iVar] = node[jPoint]->GetSolution(iVar);
  //        }
  //
  //        /*--- Set Conservative Variables ---*/
  //        numerics->SetTurbVar(Solution_i, Solution_j);
  //
  //        /*--- Retrieve flow solution for both points ---*/
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
  //          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
  //          FlowPrimVar_j[iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
  //        }
  //
  //        /*--- Set Flow Variables ---*/
  //        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
  //
  //        /*--- Set the normal vector ---*/
  //        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          Vector[iDim] = -Vector[iDim];
  //        numerics->SetNormal(Vector);
  //
  //        /*--- Add Residuals and Jacobians ---*/
  //        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
  //        LinSysRes.AddBlock(iPoint, Residual);
  //        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
  //
  //      }
  //    }
  //  }
  //
  //#else
  //
  //  int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
  //  su2double *Conserv_Var, *Flow_Var;
  //  bool compute;
  //
  //  unsigned short Buffer_Size = nVar+solver_container[FLOW_SOL]->GetnVar();
  //  su2double *Buffer_Send_U = new su2double [Buffer_Size];
  //  su2double *Buffer_Receive_U = new su2double [Buffer_Size];
  //
  //  /*--- Do the send process, by the moment we are sending each
  //   node individually, this must be changed ---*/
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
  //      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
  //      else compute = true;
  //
  //      /*--- We only send the information that belong to other boundary ---*/
  //      if ((jProcessor != rank) && compute) {
  //
  //        Conserv_Var = node[iPoint]->GetSolution();
  //        Flow_Var = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
  //
  //        for (iVar = 0; iVar < nVar; iVar++)
  //          Buffer_Send_U[iVar] = Conserv_Var[iVar];
  //
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
  //          Buffer_Send_U[nVar+iVar] = Flow_Var[iVar];
  //
  //        MPI::COMM_WORLD.Bsend(Buffer_Send_U, Buffer_Size, MPI::DOUBLE, jProcessor, iPoint);
  //
  //      }
  //    }
  //  }
  //
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
  //      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
  //      else compute = true;
  //
  //      if (compute) {
  //
  //        /*--- We only receive the information that belong to other boundary ---*/
  //        if (jProcessor != rank) {
  //          MPI::COMM_WORLD.Recv(Buffer_Receive_U, Buffer_Size, MPI::DOUBLE, jProcessor, jPoint);
  //        }
  //        else {
  //
  //          for (iVar = 0; iVar < nVar; iVar++)
  //            Buffer_Receive_U[iVar] = node[jPoint]->GetSolution(iVar);
  //
  //          for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
  //            Buffer_Send_U[nVar+iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
  //
  //        }
  //
  //        /*--- Store the solution for both points ---*/
  //        for (iVar = 0; iVar < nVar; iVar++) {
  //          Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
  //          Solution_j[iVar] = Buffer_Receive_U[iVar];
  //        }
  //
  //        /*--- Set Turbulent Variables ---*/
  //        numerics->SetTurbVar(Solution_i, Solution_j);
  //
  //        /*--- Retrieve flow solution for both points ---*/
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
  //          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
  //          FlowPrimVar_j[iVar] = Buffer_Receive_U[nVar + iVar];
  //        }
  //
  //        /*--- Set Flow Variables ---*/
  //        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
  //
  //        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          Vector[iDim] = -Vector[iDim];
  //        numerics->SetNormal(Vector);
  //
  //        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
  //        LinSysRes.AddBlock(iPoint, Residual);
  //        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
  //
  //      }
  //    }
  //  }
  //
  //  delete[] Buffer_Send_U;
  //  delete[] Buffer_Receive_U;
  //
  //#endif
  //
  //  delete[] Vector;
  //
}

void CTurbSASolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
    CNumerics *visc_numerics, CConfig *config){

  unsigned long iVertex, jVertex, iPoint, Point_Normal = 0;
  unsigned short iDim, iVar, iMarker;

  bool grid_movement = config->GetGrid_Movement();
  unsigned short nPrimVar = solver_container[FLOW_SOL]->GetnPrimVar();
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  su2double *tmp_residual = new su2double[nVar];
  
  unsigned long nDonorVertex;
  su2double weight;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

        if (geometry->node[iPoint]->GetDomain()) {
          
          nDonorVertex = GetnSlidingStates(iMarker, iVertex);
          
          /*--- Initialize Residual, this will serve to accumulate the average ---*/

          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = 0.0;

          /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/

          for (jVertex = 0; jVertex < nDonorVertex; jVertex++){

            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(iVar);
              PrimVar_j[iVar] = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, iVar, jVertex);
            }

            /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/

            weight = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

            /*--- Set primitive variables ---*/

            conv_numerics->SetPrimitive( PrimVar_i, PrimVar_j );

            /*--- Set the turbulent variable states ---*/
            Solution_i[0] = node[iPoint]->GetSolution(0);
            Solution_j[0] = GetSlidingState(iMarker, iVertex, 0, jVertex);

            conv_numerics->SetTurbVar(Solution_i, Solution_j);
            /*--- Set the normal vector ---*/

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

          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

          /*--- Set the normal vector and the coordinates ---*/

          visc_numerics->SetNormal(Normal);
          visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

          /*--- Primitive variables, and gradient ---*/

          visc_numerics->SetPrimitive(PrimVar_i, PrimVar_j);
          //          visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

          /*--- Turbulent variables and its gradients  ---*/

          visc_numerics->SetTurbVar(Solution_i, Solution_j);
          visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

          /*--- Compute and update residual ---*/

          visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

          LinSysRes.SubtractBlock(iPoint, Residual);

          /*--- Jacobian contribution for implicit integration ---*/

          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

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

void CTurbSASolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                          CConfig *config, unsigned short val_marker) {
  
  //  unsigned long iVertex, iPoint, jPoint;
  //  unsigned short iVar, iDim;
  //
  //  su2double *Vector = new su2double[nDim];
  //
  //#ifndef HAVE_MPI
  //
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();
  //
  //      if (iPoint != jPoint) {
  //
  //        /*--- Store the solution for both points ---*/
  //        for (iVar = 0; iVar < nVar; iVar++) {
  //          Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
  //          Solution_j[iVar] = node[jPoint]->GetSolution(iVar);
  //        }
  //
  //        /*--- Set Conservative Variables ---*/
  //        numerics->SetTurbVar(Solution_i, Solution_j);
  //
  //        /*--- Retrieve flow solution for both points ---*/
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
  //          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
  //          FlowPrimVar_j[iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
  //        }
  //
  //        /*--- Set Flow Variables ---*/
  //        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
  //
  //        /*--- Set the normal vector ---*/
  //        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          Vector[iDim] = -Vector[iDim];
  //        numerics->SetNormal(Vector);
  //
  //        /*--- Add Residuals and Jacobians ---*/
  //        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
  //        LinSysRes.AddBlock(iPoint, Residual);
  //        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
  //
  //      }
  //    }
  //  }
  //
  //#else
  //
  //  int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
  //  su2double *Conserv_Var, *Flow_Var;
  //  bool compute;
  //
  //  unsigned short Buffer_Size = nVar+solver_container[FLOW_SOL]->GetnVar();
  //  su2double *Buffer_Send_U = new su2double [Buffer_Size];
  //  su2double *Buffer_Receive_U = new su2double [Buffer_Size];
  //
  //  /*--- Do the send process, by the moment we are sending each
  //   node individually, this must be changed ---*/
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
  //      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
  //      else compute = true;
  //
  //      /*--- We only send the information that belong to other boundary ---*/
  //      if ((jProcessor != rank) && compute) {
  //
  //        Conserv_Var = node[iPoint]->GetSolution();
  //        Flow_Var = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
  //
  //        for (iVar = 0; iVar < nVar; iVar++)
  //          Buffer_Send_U[iVar] = Conserv_Var[iVar];
  //
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
  //          Buffer_Send_U[nVar+iVar] = Flow_Var[iVar];
  //
  //        MPI::COMM_WORLD.Bsend(Buffer_Send_U, Buffer_Size, MPI::DOUBLE, jProcessor, iPoint);
  //
  //      }
  //    }
  //  }
  //
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
  //      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
  //      else compute = true;
  //
  //      if (compute) {
  //
  //        /*--- We only receive the information that belong to other boundary ---*/
  //        if (jProcessor != rank) {
  //          MPI::COMM_WORLD.Recv(Buffer_Receive_U, Buffer_Size, MPI::DOUBLE, jProcessor, jPoint);
  //        }
  //        else {
  //
  //          for (iVar = 0; iVar < nVar; iVar++)
  //            Buffer_Receive_U[iVar] = node[jPoint]->GetSolution(iVar);
  //
  //          for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
  //            Buffer_Send_U[nVar+iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
  //
  //        }
  //
  //        /*--- Store the solution for both points ---*/
  //        for (iVar = 0; iVar < nVar; iVar++) {
  //          Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
  //          Solution_j[iVar] = Buffer_Receive_U[iVar];
  //        }
  //
  //        /*--- Set Turbulent Variables ---*/
  //        numerics->SetTurbVar(Solution_i, Solution_j);
  //
  //        /*--- Retrieve flow solution for both points ---*/
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
  //          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
  //          FlowPrimVar_j[iVar] = Buffer_Receive_U[nVar + iVar];
  //        }
  //
  //        /*--- Set Flow Variables ---*/
  //        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
  //
  //        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          Vector[iDim] = -Vector[iDim];
  //        numerics->SetNormal(Vector);
  //
  //        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
  //        LinSysRes.AddBlock(iPoint, Residual);
  //        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
  //
  //      }
  //    }
  //  }
  //
  //  delete[] Buffer_Send_U;
  //  delete[] Buffer_Receive_U;
  //
  //#endif
  //
  //  delete[] Vector;
  //
}

void CTurbSASolver::SetDES_LengthScale(CSolver **solver, CGeometry *geometry, CConfig *config){
  
  unsigned short kindHybridRANSLES = config->GetKind_HybridRANSLES();
  unsigned long iPoint = 0, jPoint = 0;
  unsigned short iDim = 0, jDim = 0, iNeigh = 0, nNeigh = 0;
  
  su2double constDES = config->GetConst_DES();
  
  su2double density = 0.0, laminarViscosity = 0.0, kinematicViscosity = 0.0,
      eddyViscosity = 0.0, kinematicViscosityTurb = 0.0, wallDistance = 0.0, lengthScale = 0.0;
  
  su2double maxDelta = 0.0, deltaAux = 0.0, distDES = 0.0, uijuij = 0.0, k2 = 0.0, r_d = 0.0, f_d = 0.0,
      deltaDDES = 0.0, deltaAuxDDES = 0.0, omega = 0.0, ln_max = 0.0, ln[3] = {0.0, 0.0, 0.0},
      aux_ln = 0.0, f_kh = 0.0;
  
  su2double nu_hat, fw_star = 0.424, cv1_3 = pow(7.1, 3.0); k2 = pow(0.41, 2.0);
  su2double cb1   = 0.1355, ct3 = 1.2, ct4   = 0.5;
  su2double sigma = 2./3., cb2 = 0.622, f_max=1.0, f_min=0.1, a1=0.15, a2=0.3;
  su2double cw1 = 0.0, Ji = 0.0, Ji_2 = 0.0, Ji_3 = 0.0, fv1 = 0.0, fv2 = 0.0, ft2 = 0.0, psi_2 = 0.0;
  su2double *coord_i = NULL, *coord_j = NULL, **primVarGrad = NULL, *vorticity = NULL, delta[3] = {0.0,0.0,0.0},
      ratioOmega[3] = {0.0, 0.0, 0.0}, vortexTiltingMeasure = 0.0;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    
    coord_i                 = geometry->node[iPoint]->GetCoord();
    nNeigh                  = geometry->node[iPoint]->GetnPoint();
    wallDistance            = geometry->node[iPoint]->GetWall_Distance();
    primVarGrad             = solver[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
    vorticity               = solver[FLOW_SOL]->node[iPoint]->GetVorticity();    
    density                 = solver[FLOW_SOL]->node[iPoint]->GetDensity();
    laminarViscosity        = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    eddyViscosity           = solver[TURB_SOL]->node[iPoint]->GetmuT();
    kinematicViscosity      = laminarViscosity/density;
    kinematicViscosityTurb  = eddyViscosity/density;
    
    uijuij = 0.0;
    for(iDim = 0; iDim < nDim; iDim++){
      for(jDim = 0; jDim < nDim; jDim++){
        uijuij += primVarGrad[1+iDim][jDim]*primVarGrad[1+iDim][jDim];
      }
    }
    uijuij = sqrt(fabs(uijuij));
    uijuij = max(uijuij,1e-10);
    
    /*--- Low Reynolds number correction term ---*/
    
    nu_hat = node[iPoint]->GetSolution()[0];
    Ji   = nu_hat/kinematicViscosity;
    Ji_2 = Ji * Ji;
    Ji_3 = Ji*Ji*Ji;
    fv1  = Ji_3/(Ji_3+cv1_3);
    fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    ft2 = ct3*exp(-ct4*Ji_2);
    cw1 = cb1/k2+(1.0+cb2)/sigma;
    
    psi_2 = (1.0 - (cb1/(cw1*k2*fw_star))*(ft2 + (1.0 - ft2)*fv2))/(fv1 * max(1.0e-10,1.0-ft2));
    psi_2 = min(100.0,psi_2);
    
    switch(kindHybridRANSLES){
      case SA_DES:
        /*--- Original Detached Eddy Simulation (DES97)
        Spalart
        1997
        ---*/
        maxDelta=0.;      
        for (iNeigh = 0;iNeigh < nNeigh; iNeigh++){
          jPoint  = geometry->node[iPoint]->GetPoint(iNeigh);
          coord_j = geometry->node[jPoint]->GetCoord();
          
          deltaAux = 0.;
          for (iDim = 0;iDim < nDim; iDim++){
            deltaAux += pow((coord_j[iDim]-coord_i[iDim]),2.);
          }
          
          maxDelta = max(maxDelta,sqrt(deltaAux));
        }
        
        distDES         = constDES * maxDelta;
        lengthScale = min(distDES,wallDistance);
                
        break;
        
      case SA_DDES:
        /*--- A New Version of Detached-eddy Simulation, Resistant to Ambiguous Grid Densities.
         Spalart et al.
         Theoretical and Computational Fluid Dynamics - 2006
         ---*/
            
        maxDelta = 0.0;      
        for (iNeigh = 0;iNeigh < nNeigh; iNeigh++){
          jPoint  = geometry->node[iPoint]->GetPoint(iNeigh);
          coord_j = geometry->node[jPoint]->GetCoord();
        
          deltaAux = 0.0;
          for (iDim = 0; iDim < nDim; iDim++){
            deltaAux += pow((coord_j[iDim]-coord_i[iDim]),2.);
          }
          
          maxDelta = max(maxDelta,sqrt(deltaAux));
        }
        
        r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
        f_d = 1.0-tanh(pow(8.0*r_d,3.0));
        
        distDES = constDES * maxDelta;
        lengthScale = wallDistance-f_d*max(0.0,(wallDistance-distDES));
        
        break;
      case SA_ZDES:
        /*--- Recent improvements in the Zonal Detached Eddy Simulation (ZDES) formulation.
         Deck
         Theoretical and Computational Fluid Dynamics - 2012
         ---*/
        
        deltaDDES = 0.0;
        for (iNeigh = 0; iNeigh < nNeigh; iNeigh++){
            jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
            coord_j = geometry->node[jPoint]->GetCoord();
            deltaAuxDDES = 0.0;
            for ( iDim = 0; iDim < nDim; iDim++){
              deltaAux       = abs(coord_j[iDim] - coord_i[iDim]);
              delta[iDim]     = max(delta[iDim], deltaAux);
              deltaAuxDDES += pow((coord_j[iDim]-coord_i[iDim]),2.);
            }
            deltaDDES = max(deltaDDES,sqrt(deltaAuxDDES));
        }
        
        omega = sqrt(vorticity[0]*vorticity[0] + 
                     vorticity[1]*vorticity[1] +
                     vorticity[2]*vorticity[2]);
        
        for (iDim = 0; iDim < 3; iDim++){
          ratioOmega[iDim] = vorticity[iDim]/omega;
        }
  
        maxDelta = sqrt(pow(ratioOmega[0],2.0)*delta[1]*delta[2] +
                        pow(ratioOmega[1],2.0)*delta[0]*delta[2] +
                        pow(ratioOmega[2],2.0)*delta[0]*delta[1]);
            
        r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
        f_d = 1.0-tanh(pow(8.0*r_d,3.0));
        
        if (f_d < 0.99){
          maxDelta = deltaDDES;
        }
        
        distDES = constDES * maxDelta;
        lengthScale = wallDistance-f_d*max(0.0,(wallDistance-distDES));
        
        break;
        
      case SA_EDDES:
        
        /*--- An Enhanced Version of DES with Rapid Transition from RANS to LES in Separated Flows.
         Shur et al.
         Flow Turbulence Combust - 2015
         ---*/
        
        vortexTiltingMeasure = node[iPoint]->GetVortex_Tilting();
        
        omega = sqrt(vorticity[0]*vorticity[0] + 
                     vorticity[1]*vorticity[1] +
                     vorticity[2]*vorticity[2]);
        
        for (iDim = 0; iDim < 3; iDim++){
          ratioOmega[iDim] = vorticity[iDim]/omega;
        }
        
        ln_max = 0.0;
        deltaDDES = 0.0;
        for (iNeigh = 0;iNeigh < nNeigh; iNeigh++){
          jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
          coord_j = geometry->node[jPoint]->GetCoord();
          deltaAuxDDES = 0.0;
          for (iDim = 0; iDim < nDim; iDim++){
            delta[iDim] = fabs(coord_j[iDim] - coord_i[iDim]);            
            deltaAuxDDES += pow((coord_j[iDim]-coord_i[iDim]),2.);
          }
          deltaDDES=max(deltaDDES,sqrt(deltaAuxDDES));
          ln[0] = delta[1]*ratioOmega[2] - delta[2]*ratioOmega[1];
          ln[1] = delta[2]*ratioOmega[0] - delta[0]*ratioOmega[2];
          ln[2] = delta[0]*ratioOmega[1] - delta[1]*ratioOmega[0];
          aux_ln = sqrt(ln[0]*ln[0] + ln[1]*ln[1] + ln[2]*ln[2]);
          ln_max = max(ln_max,aux_ln);
          vortexTiltingMeasure += node[jPoint]->GetVortex_Tilting();
        }

        vortexTiltingMeasure = (vortexTiltingMeasure/fabs(nNeigh + 1.0));
        
        f_kh = max(f_min, min(f_max, f_min + ((f_max - f_min)/(a2 - a1)) * (vortexTiltingMeasure - a1)));
        
        r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
        f_d = 1.0-tanh(pow(8.0*r_d,3.0));

        maxDelta = (ln_max/sqrt(3.0)) * f_kh;
        if (f_d < 0.999){
          maxDelta = deltaDDES;
        }
        
        distDES = constDES * maxDelta;
        lengthScale=wallDistance-f_d*max(0.0,(wallDistance-distDES));
        
        break;
        
    }
    
    node[iPoint]->SetDES_LengthScale(lengthScale);
  
  }
}

void CTurbSASolver::SetInlet(CConfig* config) {

  /*-- Allocation has to happen in derived classes (not CTurbSolver),
   * due to arbitrary number of turbulence variables ---*/
  Inlet_TurbVars = new su2double**[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    bool isCustomizable = config->GetMarker_All_PyCustom(iMarker);
    bool isInlet = (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW);
    if (isCustomizable && isInlet) {
      Inlet_TurbVars[iMarker] = new su2double*[nVertex[iMarker]];
      for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
        Inlet_TurbVars[iMarker][iVertex] = new su2double[nVar];
        /*--- The default values for these custom BCs are initialized with
         * the freestream values to avoid nonphysical values at unspecified
         * points ---*/
        Inlet_TurbVars[iMarker][iVertex][0] = nu_tilde_Inf;
      }
    } else {
      Inlet_TurbVars[iMarker] = NULL;
    }
  }

}

CTurbSSTSolver::CTurbSSTSolver(void) : CTurbSolver() {
  
  /*--- Array initialization ---*/
  constants = NULL;
  Inlet_TurbVars = NULL;
  
}

CTurbSSTSolver::CTurbSSTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver(geometry, config) {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;
  
  /*--- Array initialization ---*/
  
  constants = NULL;
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Dimension of the problem --> dependent on the turbulence model. ---*/
  
  nVar = 2;
  nPrimVar = 2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  
  /*--- Single grid simulation ---*/
  
  if (iMesh == MESH_0) {
    
    /*--- Define some auxiliary vector related with the residual ---*/
    
    Residual = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]  = 0.0;
    Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_i = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]  = 0.0;
    Residual_j = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]  = 0.0;
    Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
    
    /*--- Define some structures for locating max residuals ---*/
    
    Point_Max = new unsigned long[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
    Point_Max_Coord = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
    }
    
    /*--- Define some auxiliary vector related with the solution ---*/
    
    Solution = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
    
    /*--- Define some auxiliary vector related with the geometry ---*/
    
    Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
    
    /*--- Define some auxiliary vector related with the flow solution ---*/
    
    FlowPrimVar_i = new su2double [nDim+7]; FlowPrimVar_j = new su2double [nDim+7];
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    /*--- Initialization of the structure of the whole Jacobian ---*/
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SST model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }
  
  /*--- Computation of gradients by least squares ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
    Smatrix[iDim] = new su2double [nDim];
    /*--- c vector := transpose(WA)*(Wb) ---*/
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
    Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Initialize value for model constants ---*/
  constants = new su2double[10];
  constants[0] = 0.85;   //sigma_k1
  constants[1] = 1.0;    //sigma_k2
  constants[2] = 0.5;    //sigma_om1
  constants[3] = 0.856;  //sigma_om2
  constants[4] = 0.075;  //beta_1
  constants[5] = 0.0828; //beta_2
  constants[6] = 0.09;   //betaStar
  constants[7] = 0.31;   //a1
  constants[8] = constants[4]/constants[6] - constants[2]*0.41*0.41/sqrt(constants[6]);  //alfa_1
  constants[9] = constants[5]/constants[6] - constants[3]*0.41*0.41/sqrt(constants[6]);  //alfa_2
  
  /*--- Initialize lower and upper limits---*/
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  
  lowerlimit[0] = 1.0e-10;
  upperlimit[0] = 1.0e10;
  
  lowerlimit[1] = 1.0e-4;
  upperlimit[1] = 1.0e15;
  
  /*--- Far-field flow state quantities and initialization. ---*/
  su2double rhoInf, *VelInf, muLamInf, Intensity, viscRatio, muT_Inf;

  rhoInf    = config->GetDensity_FreeStreamND();
  VelInf    = config->GetVelocity_FreeStreamND();
  muLamInf  = config->GetViscosity_FreeStreamND();
  Intensity = config->GetTurbulenceIntensity_FreeStream();
  viscRatio = config->GetTurb2LamViscRatio_FreeStream();

  su2double VelMag = 0;
  for (iDim = 0; iDim < nDim; iDim++)
    VelMag += VelInf[iDim]*VelInf[iDim];
  VelMag = sqrt(VelMag);

  kine_Inf  = 3.0/2.0*(VelMag*VelMag*Intensity*Intensity);
  omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);

  /*--- Eddy viscosity, initialized without stress limiter at the infinity ---*/
  muT_Inf = rhoInf*kine_Inf/omega_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CTurbSSTVariable(kine_Inf, omega_Inf, muT_Inf, nDim, nVar, constants, config);

  /*--- MPI solution ---*/

//TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart
  Set_MPI_Solution(geometry, config);
  Set_MPI_Solution(geometry, config);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  unsigned long iMarker;

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

  /*--- Set up inlet profiles, if necessary ---*/

  SetInlet(config);

}

CTurbSSTSolver::~CTurbSSTSolver(void) {
  
  if (constants != NULL) delete [] constants;
  
  unsigned long iMarker, iVertex;
  unsigned short iVar;

  if ( SlidingState != NULL ) {
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

}

void CTurbSSTSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint;

  unsigned long ExtIter = config->GetExtIter();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter_flow     = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_turb     = ((config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the residual vector ---*/
    
    LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  Jacobian.SetValZero();

  /*--- Upwind second order reconstruction ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

  if (limiter_turb) SetSolution_Limiter(geometry, config);
  
  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);

}

void CTurbSSTSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  su2double rho = 0.0, mu = 0.0, dist, omega, kine, strMag, F2, muT, zeta;
  su2double a1 = constants[7];
  unsigned long iPoint;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Compute mean flow and turbulence gradients ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
//    solver_container[FLOW_SOL]->SetPrimitive_Gradient_GG(geometry, config);
    SetSolution_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
//    solver_container[FLOW_SOL]->SetPrimitive_Gradient_LS(geometry, config);
    SetSolution_Gradient_LS(geometry, config);
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Compute blending functions and cross diffusion ---*/
    
    if (compressible) {
      rho  = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      mu   = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    }
    if (incompressible) {
      rho  = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      mu   = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    }
    
    dist = geometry->node[iPoint]->GetWall_Distance();
    
    strMag = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();

    node[iPoint]->SetBlendingFunc(mu, dist, rho);
    
    F2 = node[iPoint]->GetF2blending();
    
    /*--- Compute the eddy viscosity ---*/
    
    kine  = node[iPoint]->GetSolution(0);
    omega = node[iPoint]->GetSolution(1);
    zeta = min(1.0/omega, a1/(strMag*F2));
    muT = min(max(rho*kine*zeta,0.0),1.0);
    node[iPoint]->SetmuT(muT);
    
  }
  
}

void CTurbSSTSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {
  
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);
    
    /*--- Gradient of the primitive and conservative variables ---*/
    
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
    
    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    
    numerics->SetTurbVar(node[iPoint]->GetSolution(), NULL);
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);
    
    /*--- Set volume ---*/
    
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- Set distance to the surface ---*/
    
    numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
    
    /*--- Menter's first blending function ---*/
    
    numerics->SetF1blending(node[iPoint]->GetF1blending(),0.0);
    
    /*--- Menter's second blending function ---*/
    
    numerics->SetF2blending(node[iPoint]->GetF2blending(),0.0);
    
    /*--- Set vorticity and strain rate magnitude ---*/
    
    numerics->SetVorticity(solver_container[FLOW_SOL]->node[iPoint]->GetVorticity(), NULL);
    
    numerics->SetStrainMag(solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag(), 0.0);
    
    /*--- Cross diffusion ---*/
    
    numerics->SetCrossDiff(node[iPoint]->GetCrossDiff(),0.0);
    
    /*--- Compute the source term ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);
    
    /*--- Subtract residual and the Jacobian ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    
  }
  
}

void CTurbSSTSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {
  
}

void CTurbSSTSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;
  su2double distance, density = 0.0, laminar_viscosity = 0.0, beta_1;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- distance to closest neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      distance = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
        (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      distance = sqrt(distance);
      
      /*--- Set wall values ---*/
      if (compressible) {
        density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
        laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      }
      if (incompressible) {
        density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
        laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      }
      
      beta_1 = constants[4];
      
      Solution[0] = 0.0;
      Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);
      
      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution);
      node[iPoint]->SetSolution(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      
    }
  }
  
}

void CTurbSSTSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                        unsigned short val_marker) {
  
  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;
  su2double distance, density = 0.0, laminar_viscosity = 0.0, beta_1;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- distance to closest neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      distance = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
        (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      distance = sqrt(distance);
      
      /*--- Set wall values ---*/
      if (compressible) {
        density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
        laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      }
      if (incompressible) {
        density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
        laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      }
      
      beta_1 = constants[4];
      
      Solution[0] = 0.0;
      Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);
      
      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution);
      node[iPoint]->SetSolution(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      
    }
  }
  
}

void CTurbSSTSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  su2double *Normal, *V_infty, *V_domain;
  unsigned short iVar, iDim;
  
  bool grid_movement = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set turbulent variable at the wall, and at infinity ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
      Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      
      Solution_j[0] = kine_Inf;
      Solution_j[1] = omega_Inf;
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute residuals and Jacobians ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  delete [] Normal;
  
}

void CTurbSSTSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker) {
  
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint;
  su2double *V_inlet, *V_domain, *Normal;
  bool isCustomizable = config->GetMarker_All_PyCustom(val_marker);
  
  Normal = new su2double[nDim];
  
  bool grid_movement  = config->GetGrid_Movement();
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the inlet ---*/
      
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Set the turbulent variable states. Use free-stream SST
       values for the turbulent state at the inflow. ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
      Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      
      if (isCustomizable) {
        /*--- Only use the array version of the inlet BC when necessary ---*/
        Solution_j[0] = Inlet_TurbVars[val_marker][iVertex][0];
        Solution_j[1] = Inlet_TurbVars[val_marker][iVertex][1];
      } else {
        Solution_j[0] = kine_Inf;
        Solution_j[1] = omega_Inf;
      }
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_inlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//     visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- Menter's first blending function ---*/
//
//      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }
      
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

void CTurbSSTSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the outlet ---*/
      
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TurbVar_internal,
       Solution_j --> TurbVar_outlet ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- Menter's first blending function ---*/
//
//      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}


void CTurbSSTSolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker) {

  unsigned short iVar, iSpan, iDim;
  unsigned long  oldVertex, iPoint, Point_Normal;
  long iVertex;
  su2double *V_inlet, *V_domain, *Normal;
  su2double extAverageKine, extAverageOmega;
  unsigned short nSpanWiseSections = config->GetnSpanWiseSections();

  Normal = new su2double[nDim];

  bool grid_movement  = config->GetGrid_Movement();

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){
    extAverageKine = solver_container[FLOW_SOL]->GetExtAverageKine(val_marker, iSpan);
    extAverageOmega = solver_container[FLOW_SOL]->GetExtAverageOmega(val_marker, iSpan);


    /*--- Loop over all the vertices on this boundary marker ---*/

    for (iVertex = 0; iVertex < geometry->nVertexSpan[val_marker][iSpan]; iVertex++) {

      /*--- find the node related to the vertex ---*/
      iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);

      Solution_j[0]= extAverageKine;
      Solution_j[1]= extAverageOmega;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
            geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());

      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}

void CTurbSSTSolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker) {

  unsigned short iVar, iSpan, iDim;
  unsigned long  oldVertex, iPoint, Point_Normal;
  long iVertex;
  su2double *V_inlet, *V_domain, *Normal;
  unsigned short nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Quantities for computing the  kine and omega to impose at the inlet boundary. ---*/
  su2double rho, pressure, *Vel, VelMag, muLam, Intensity, viscRatio, kine_b, omega_b, kine;
  CFluidModel *FluidModel;

  FluidModel = solver_container[FLOW_SOL]->GetFluidModel();
  Intensity = config->GetTurbulenceIntensity_FreeStream();
  viscRatio = config->GetTurb2LamViscRatio_FreeStream();

  Normal = new su2double[nDim];
  Vel = new su2double[nDim];

  bool grid_movement  = config->GetGrid_Movement();

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);


  for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){

    /*--- Compute the inflow kine and omega using the span wise averge quntities---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Vel[iDim] = solver_container[FLOW_SOL]->GetAverageTurboVelocity(val_marker, iSpan)[iDim];

    rho       = solver_container[FLOW_SOL]->GetAverageDensity(val_marker, iSpan);
    pressure  = solver_container[FLOW_SOL]->GetAveragePressure(val_marker, iSpan);
    kine      = solver_container[FLOW_SOL]->GetAverageKine(val_marker, iSpan);

    FluidModel->SetTDState_Prho(pressure, rho);
    muLam = FluidModel->GetLaminarViscosity();

    VelMag = 0;
    for (iDim = 0; iDim < nDim; iDim++)
      VelMag += Vel[iDim]*Vel[iDim];
    VelMag = sqrt(VelMag);

    kine_b  = 3.0/2.0*(VelMag*VelMag*Intensity*Intensity);
    omega_b = rho*kine/(muLam*viscRatio);

    /*--- Loop over all the vertices on this boundary marker ---*/
    for (iVertex = 0; iVertex < geometry->nVertexSpan[val_marker][iSpan]; iVertex++) {

      /*--- find the node related to the vertex ---*/
      iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);

      /*--- Set the turbulent variable states. Use average span-wise values
             values for the turbulent state at the inflow. ---*/

      Solution_j[0]= kine_b;
      Solution_j[1]= omega_b;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
            geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());

      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  delete[] Vel;

}


void CTurbSSTSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
    CNumerics *visc_numerics, CConfig *config){

  unsigned long iVertex, jVertex, iPoint, Point_Normal = 0;
  unsigned short iDim, iVar, iMarker;

  bool grid_movement = config->GetGrid_Movement();
  unsigned short nPrimVar = solver_container[FLOW_SOL]->GetnPrimVar();
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  su2double *tmp_residual = new su2double[nVar];
  
  unsigned long nDonorVertex;
  su2double weight;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

        if (geometry->node[iPoint]->GetDomain()) {

          nDonorVertex = GetnSlidingStates(iMarker, iVertex);

          /*--- Initialize Residual, this will serve to accumulate the average ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = 0.0;

          /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/
          
          for (jVertex = 0; jVertex < nDonorVertex; jVertex++){
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(iVar);
              PrimVar_j[iVar] = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, iVar, jVertex);
            }

            /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/
            
            weight = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

            /*--- Set primitive variables ---*/

            conv_numerics->SetPrimitive( PrimVar_i, PrimVar_j );

            /*--- Set the turbulent variable states ---*/
            Solution_i[0] = node[iPoint]->GetSolution(0);
            Solution_i[1] = node[iPoint]->GetSolution(1);

            Solution_j[0] = GetSlidingState(iMarker, iVertex, 0, jVertex);
            Solution_j[1] = GetSlidingState(iMarker, iVertex, 1, jVertex);

            conv_numerics->SetTurbVar(Solution_i, Solution_j);
    
            /*--- Set the normal vector ---*/

            conv_numerics->SetNormal(Normal);

            if (grid_movement)
              conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

            conv_numerics->ComputeResidual(tmp_residual, Jacobian_i, Jacobian_j, config);

            /*--- Accumulate the residuals to compute the average ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
              Residual[iVar] += weight*tmp_residual[iVar];
          }
          
          /*--- Add Residuals and Jacobians ---*/

          LinSysRes.AddBlock(iPoint, Residual);

          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

          /*--- Set the normal vector and the coordinates ---*/

          visc_numerics->SetNormal(Normal);
          visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

          /*--- Primitive variables, and gradient ---*/

          visc_numerics->SetPrimitive(PrimVar_i, PrimVar_j);
          //          visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

          /*--- Turbulent variables and its gradients  ---*/

          visc_numerics->SetTurbVar(Solution_i, Solution_j);
          visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

          /*--- Compute and update residual ---*/

          visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

          LinSysRes.SubtractBlock(iPoint, Residual);

          /*--- Jacobian contribution for implicit integration ---*/

          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

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

su2double* CTurbSSTSolver::GetConstants() {
  return constants;
}


void CTurbSSTSolver::SetInlet(CConfig* config) {

  /*-- Allocation has to happen in derived classes (not CTurbSolver),
   * due to arbitrary number of turbulence variables ---*/
  Inlet_TurbVars = new su2double**[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    bool isCustomizable = config->GetMarker_All_PyCustom(iMarker);
    bool isInlet = (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW);
    if (isCustomizable && isInlet) {
      Inlet_TurbVars[iMarker] = new su2double*[nVertex[iMarker]];
      for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
        Inlet_TurbVars[iMarker][iVertex] = new su2double[nVar];
        /*--- The default values for these custom BCs are initialized with
         * the freestream values to avoid nonphysical values at unspecified
         * points ---*/
        Inlet_TurbVars[iMarker][iVertex][0] = kine_Inf;
        Inlet_TurbVars[iMarker][iVertex][1] = omega_Inf;
      }
    } else {
      Inlet_TurbVars[iMarker] = NULL;
    }
  }

}
