/*!
 * \file element_structure.inl
 * \brief In-Line subroutines of the <i>element_structure.hpp</i> file.
 * \author R. Sanchez
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
 
#pragma once

inline void CElement::ComputeGrad_Linear(void) { }

inline void CElement::ComputeGrad_NonLinear(void) { }

inline unsigned short CElement::GetnNodes(void) { return nNodes;}

inline unsigned short CElement::GetnGaussPoints(void) { return nGaussPoints;}

inline void CElement::SetRef_Coord(su2double val_CoordRef, unsigned short iNode, unsigned short iDim) { RefCoord[iNode][iDim] = val_CoordRef;}

inline void CElement::SetCurr_Coord(su2double val_CoordCurr, unsigned short iNode, unsigned short iDim) { CurrentCoord[iNode][iDim] = val_CoordCurr;}

inline void CElement::SetElement_Pressure(su2double val_ElPressure) {el_Pressure = val_ElPressure;}

inline su2double CElement::GetRef_Coord(unsigned short iNode, unsigned short iDim) { return RefCoord[iNode][iDim];}

inline su2double CElement::GetCurr_Coord(unsigned short iNode, unsigned short iDim) { return CurrentCoord[iNode][iDim];}

inline su2double CElement::GetWeight(unsigned short iGauss) { return GaussWeight[iGauss];}

inline su2double CElement::GetJ_X(unsigned short iGauss) { return GaussPoint[iGauss]->GetJ_X();}

inline su2double CElement::GetJ_x(unsigned short iGauss) { return GaussPoint[iGauss]->GetJ_x();}

inline su2double CElement::GetJ_X_P(unsigned short iGauss) { return GaussPointP[iGauss]->GetJ_X();}

inline su2double CElement::GetJ_x_P(unsigned short iGauss) { return GaussPointP[iGauss]->GetJ_x();}

inline su2double CElement::GetElement_Pressure(void) { return el_Pressure;}

inline su2double CElement::Get_Mab(unsigned short nodeA, unsigned short nodeB) { return Mab[nodeA][nodeB]; }

inline su2double *CElement::Get_Kab(unsigned short nodeA, unsigned short nodeB) { return Kab[nodeA][nodeB];}

inline su2double *CElement::Get_Kk_ab(unsigned short nodeA, unsigned short nodeB) { return Kk_ab[nodeA][nodeB];}

inline su2double *CElement::Get_Kt_a(unsigned short nodeA) { return Kt_a[nodeA];}

inline su2double *CElement::Get_FDL_a(unsigned short nodeA) { return FDL_a[nodeA];}

inline su2double CElement::Get_Ks_ab(unsigned short nodeA, unsigned short nodeB) { return Ks_ab[nodeA][nodeB]; }

inline void CElement::Add_Mab(su2double val_Mab, unsigned short nodeA, unsigned short nodeB) { Mab[nodeA][nodeB] += val_Mab; }

inline void CElement::Add_Ks_ab(su2double val_Ks_ab, unsigned short nodeA, unsigned short nodeB) { Ks_ab[nodeA][nodeB] += val_Ks_ab; }

inline void CElement::Add_NodalStress(su2double val_Stress, unsigned short iNode, unsigned short iVar) { NodalStress[iNode][iVar] += val_Stress; }

inline su2double CElement::GetNi(unsigned short iNode, unsigned short iGauss) { return GaussPoint[iGauss]->GetNi(iNode);}

inline su2double CElement::GetGradNi_X(unsigned short iNode, unsigned short iGauss, unsigned short iDim) { return GaussPoint[iGauss]->GetGradNi_Xj(iNode,iDim);}

inline su2double CElement::GetGradNi_x(unsigned short iNode, unsigned short iGauss, unsigned short iDim) { return GaussPoint[iGauss]->GetGradNi_xj(iNode,iDim);}

inline su2double CElement::GetGradNi_X_P(unsigned short iNode, unsigned short iGaussP, unsigned short iDim) { return GaussPointP[iGaussP]->GetGradNi_Xj(iNode,iDim);}

inline su2double CElement::GetGradNi_x_P(unsigned short iNode, unsigned short iGaussP, unsigned short iDim) { return GaussPointP[iGaussP]->GetGradNi_xj(iNode,iDim);}

inline su2double CElement::GetNi_Extrap(unsigned short iNode, unsigned short iGauss) { return NodalExtrap[iNode][iGauss]; }

inline su2double CElement::Get_NodalStress(unsigned short iNode, unsigned short iVar) { return NodalStress[iNode][iVar]; }

inline su2double CElement::GetWeight_P(unsigned short iGaussP) { return GaussWeightP[iGaussP];}

inline su2double CElement::ComputeArea(void) { return 0.0;}

inline su2double CElement::ComputeVolume(void) { return 0.0;}

inline su2double CElement::ComputeCurrentArea(void) { return 0.0;}

inline su2double CElement::ComputeCurrentVolume(void) { return 0.0;}

inline unsigned short CElement::GetnGaussPointsP(void) { return nGaussPointsP;}

inline void CElement::Set_iDe(unsigned short val_iDe) { iDe = val_iDe;}

inline unsigned short CElement::Get_iDe(void) { return iDe;}

inline unsigned long CElement::Get_iDV(void) { return iDV;}

inline unsigned long CElement::Get_iProp(void) { return iProp;}

inline void CElement::ComputeGrad_Pressure(void){ }

inline void CQUAD4::ComputeGrad_Pressure(void) { }

inline void CHEXA8::ComputeGrad_Pressure(void) { }

