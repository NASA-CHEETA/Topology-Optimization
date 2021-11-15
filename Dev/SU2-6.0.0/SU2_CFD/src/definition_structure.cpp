/*!
 * \file definition_structure.cpp
 * \brief Main subroutines used by SU2_CFD
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

#include "../include/definition_structure.hpp"


void Partition_Analysis(CGeometry *geometry, CConfig *config) {
  
  /*--- This routine does a quick and dirty output of the total
   vertices, ghost vertices, total elements, ghost elements, etc.,
   so that we can analyze the partition quality. ---*/
  
  unsigned short nMarker = config->GetnMarker_All();
  unsigned short iMarker, iNodes, MarkerS, MarkerR;
  unsigned long iElem, iPoint, nVertexS, nVertexR;
  unsigned long nNeighbors = 0, nSendTotal = 0, nRecvTotal = 0;
  unsigned long nPointTotal=0, nPointGhost=0, nElemTotal=0;
  unsigned long nElemHalo=0, nEdge=0, nElemBound=0;
  int iRank;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  SU2_MPI::Comm_rank(MPI_COMM_WORLD, &rank);
  SU2_MPI::Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  nPointTotal = geometry->GetnPoint();
  nPointGhost = geometry->GetnPoint() - geometry->GetnPointDomain();
  nElemTotal  = geometry->GetnElem();
  nEdge       = geometry->GetnEdge();
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    nElemBound  += geometry->GetnElem_Bound(iMarker);
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      MarkerS = iMarker;  MarkerR = iMarker+1;
      nVertexS = geometry->nVertex[MarkerS];
      nVertexR = geometry->nVertex[MarkerR];
      nNeighbors++;
      nSendTotal += nVertexS;
      nRecvTotal += nVertexR;
    }
  }
  
  bool *isHalo = new bool[geometry->GetnElem()];
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    isHalo[iElem] = false;
    for (iNodes = 0; iNodes < geometry->elem[iElem]->GetnNodes(); iNodes++) {
      iPoint = geometry->elem[iElem]->GetNode(iNodes);
      if (!geometry->node[iPoint]->GetDomain()) isHalo[iElem] = true;
    }
  }
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if (isHalo[iElem]) nElemHalo++;
  }
  
  unsigned long *row_ptr = NULL, nnz;
  unsigned short *nNeigh = NULL;
  vector<unsigned long>::iterator it;
  vector<unsigned long> vneighs;
  
  /*--- Don't delete *row_ptr, *col_ind because they are
   asigned to the Jacobian structure. ---*/
  
  /*--- Compute the number of neighbors ---*/
  
  nNeigh = new unsigned short [geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    // +1 -> to include diagonal element
    nNeigh[iPoint] = (geometry->node[iPoint]->GetnPoint()+1);
  }
  
  /*--- Create row_ptr structure, using the number of neighbors ---*/
  
  row_ptr = new unsigned long [geometry->GetnPoint()+1];
  row_ptr[0] = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    row_ptr[iPoint+1] = row_ptr[iPoint] + nNeigh[iPoint];
  nnz = row_ptr[geometry->GetnPoint()];
  
  delete [] row_ptr;
  delete [] nNeigh;
  
  /*--- Now put this info into a CSV file for processing ---*/
  
  char cstr[200];
  ofstream Profile_File;
  strcpy (cstr, "partitioning.csv");
  Profile_File.precision(15);
  
  if (rank == MASTER_NODE) {
    /*--- Prepare and open the file ---*/
    Profile_File.open(cstr, ios::out);
    /*--- Create the CSV header ---*/
    Profile_File << "\"Rank\", \"nNeighbors\", \"nPointTotal\", \"nEdge\", \"nPointGhost\", \"nSendTotal\", \"nRecvTotal\", \"nElemTotal\", \"nElemBoundary\", \"nElemHalo\", \"nnz\"" << endl;
    Profile_File.close();
  }
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- Loop through the map and write the results to the file ---*/
  
  for (iRank = 0; iRank < size; iRank++) {
    if (rank == iRank) {
      Profile_File.open(cstr, ios::out | ios::app);
      Profile_File << rank << ", " << nNeighbors << ", " << nPointTotal << ", " << nEdge << "," << nPointGhost << ", " << nSendTotal << ", " << nRecvTotal << ", " << nElemTotal << "," << nElemBound << ", " << nElemHalo << ", " << nnz << endl;
      Profile_File.close();
    }
#ifdef HAVE_MPI
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }
  
  delete [] isHalo;
  
}
