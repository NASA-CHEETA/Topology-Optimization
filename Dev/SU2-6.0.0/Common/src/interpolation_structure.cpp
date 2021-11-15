/*!
 * \file interpolation_structure.cpp
 * \brief Main subroutines used by SU2_FSI
 * \author H. Kline
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

#include "../include/interpolation_structure.hpp"

CInterpolator::CInterpolator(void) {
  
  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  nZone = 0;
  Geometry = NULL;

  donor_geometry  = NULL;
  target_geometry = NULL;

  donorZone  = 0;
  targetZone = 0;

  Buffer_Receive_nVertex_Donor     = NULL;
  Buffer_Receive_nFace_Donor       = NULL;
  Buffer_Receive_nFaceNodes_Donor  = NULL;
  Buffer_Send_nVertex_Donor        = NULL;
  Buffer_Send_nFace_Donor          = NULL;
  Buffer_Send_nFaceNodes_Donor     = NULL;
  Buffer_Receive_GlobalPoint       = NULL;
  Buffer_Send_GlobalPoint          = NULL;
  Buffer_Send_FaceIndex            = NULL;
  Buffer_Receive_FaceIndex         = NULL;
  Buffer_Send_FaceNodes            = NULL;
  Buffer_Receive_FaceNodes         = NULL;
  Buffer_Send_FaceProc             = NULL;
  Buffer_Receive_FaceProc          = NULL;

  Buffer_Send_Coord                = NULL;
  Buffer_Send_Normal               = NULL;
  Buffer_Receive_Coord             = NULL;
  Buffer_Receive_Normal            = NULL;
  
  Receive_GlobalPoint              = NULL;
  Buffer_Receive_nLinkedNodes      = NULL;
  Buffer_Receive_LinkedNodes       = NULL;
  Buffer_Receive_StartLinkedNodes  = NULL;
  Buffer_Receive_Proc              = NULL;   

}

CInterpolator::~CInterpolator(void) {

  //if (Buffer_Receive_nVertex_Donor!= NULL) delete[] Buffer_Receive_nVertex_Donor;
}


CInterpolator::CInterpolator(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone) {
  
  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  /* Store pointers*/
  Geometry = geometry_container;

  donorZone  = iZone;
  targetZone = jZone;

  donor_geometry  = geometry_container[donorZone][MESH_0];
  target_geometry = geometry_container[targetZone][MESH_0];

  /*--- Initialize transfer coefficients between the zones ---*/
    /* Since this is a virtual function, call it in the child class constructor  */
  //Set_TransferCoeff(targetZone,donorZone,config);
  /*--- Initialize transfer coefficients between the zones ---*/
  //Set_TransferCoeff(Zones,config);

  //Buffer_Receive_nVertex_Donor = NULL;

}

inline void CInterpolator::Set_TransferCoeff(CConfig **config) { }

void CInterpolator::Determine_ArraySize(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim) {
  unsigned long nLocalVertex_Donor = 0, nLocalFaceNodes_Donor=0, nLocalFace_Donor=0;
  unsigned long iVertex, iPointDonor = 0;
  /* Only needed if face data is also collected */
  unsigned long inode;
  unsigned long donor_elem, jElem, jPoint;
  unsigned short iDonor;
  unsigned int nFaces=0, iFace, nNodes=0;
  bool face_on_marker = true;

  for (iVertex = 0; iVertex < nVertexDonor; iVertex++) {
    iPointDonor = donor_geometry->vertex[markDonor][iVertex]->GetNode();
    if (donor_geometry->node[iPointDonor]->GetDomain()) {
      nLocalVertex_Donor++;
      if (faces) {
        /*--- On Donor geometry also communicate face info ---*/
        if (nDim==3) {
          for (jElem=0; jElem<donor_geometry->node[iPointDonor]->GetnElem(); jElem++) {
            donor_elem = donor_geometry->node[iPointDonor]->GetElem(jElem);
            nFaces = donor_geometry->elem[donor_elem]->GetnFaces();
            for (iFace=0; iFace<nFaces; iFace++) {
              face_on_marker=true;
              nNodes = donor_geometry->elem[donor_elem]->GetnNodesFace(iFace);
              for (iDonor=0; iDonor<nNodes; iDonor++) {
                /*--- Local index of the node on face --*/
                inode = donor_geometry->elem[donor_elem]->GetFaces(iFace, iDonor);
                jPoint = donor_geometry->elem[donor_elem]->GetNode(inode);
                face_on_marker = (face_on_marker && (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
              }
              if (face_on_marker ) {
                nLocalFace_Donor++;
                nLocalFaceNodes_Donor+=nNodes;
              }
            }
          }
        }
        else {
          /*--- in 2D we use the edges ---*/
          nNodes=2;
          nFaces = donor_geometry->node[iPointDonor]->GetnPoint();
          for (iFace=0; iFace<nFaces; iFace++) {
            face_on_marker=true;
            for (iDonor=0; iDonor<nNodes; iDonor++) {
              inode = donor_geometry->node[iPointDonor]->GetEdge(iFace);
              jPoint = donor_geometry->edge[inode]->GetNode(iDonor);
              face_on_marker = (face_on_marker && (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
            }
            if (face_on_marker ) {
              nLocalFace_Donor++;
              nLocalFaceNodes_Donor+=nNodes;
            }
          }
        }
      }
    }
  }

  Buffer_Send_nVertex_Donor[0] = nLocalVertex_Donor;
  if (faces) {
    Buffer_Send_nFace_Donor[0] = nLocalFace_Donor;
    Buffer_Send_nFaceNodes_Donor[0] = nLocalFaceNodes_Donor;
  }

  /*--- Send Interface vertex information --*/
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalVertex_Donor, &MaxLocalVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  if (faces) {
    SU2_MPI::Allreduce(&nLocalFace_Donor, &nGlobalFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor, &MaxFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &nGlobalFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &MaxFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFace_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MaxFace_Donor++;
  }
#else
  MaxLocalVertex_Donor    = nLocalVertex_Donor;
  Buffer_Receive_nVertex_Donor[0] = Buffer_Send_nVertex_Donor[0];
  if (faces) {
    nGlobalFace_Donor       = nLocalFace_Donor;
    nGlobalFaceNodes_Donor  = nLocalFaceNodes_Donor;
    MaxFaceNodes_Donor      = nLocalFaceNodes_Donor;
    MaxFace_Donor           = nLocalFace_Donor+1;
    Buffer_Receive_nFace_Donor[0] = Buffer_Send_nFace_Donor[0];
    Buffer_Receive_nFaceNodes_Donor[0] = Buffer_Send_nFaceNodes_Donor[0];
  }
#endif

}

void CInterpolator::Collect_VertexInfo(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim) {
  unsigned long nLocalVertex_Donor = 0;
  unsigned long iVertex, iPointDonor = 0, iVertexDonor, nBuffer_Coord, nBuffer_Point;
  unsigned short iDim;
  /* Only needed if face data is also collected */
  su2double  *Normal;

  for (iVertex = 0; iVertex < MaxLocalVertex_Donor; iVertex++) {
    Buffer_Send_GlobalPoint[iVertex] = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
      if (faces)
        Buffer_Send_Normal[iVertex*nDim+iDim] = 0.0;
    }
  }

  /*--- Copy coordinates and point to the auxiliar vector --*/
  nLocalVertex_Donor = 0;

  for (iVertexDonor = 0; iVertexDonor < nVertexDonor; iVertexDonor++) {
    iPointDonor = donor_geometry->vertex[markDonor][iVertexDonor]->GetNode();
    if (donor_geometry->node[iPointDonor]->GetDomain()) {
      Buffer_Send_GlobalPoint[nLocalVertex_Donor] = donor_geometry->node[iPointDonor]->GetGlobalIndex();
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[nLocalVertex_Donor*nDim+iDim] = donor_geometry->node[iPointDonor]->GetCoord(iDim);

      if (faces) {
        Normal =  donor_geometry->vertex[markDonor][iVertexDonor]->GetNormal();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Normal[nLocalVertex_Donor*nDim+iDim] = Normal[iDim];
      }
      nLocalVertex_Donor++;
    }
  }
  nBuffer_Coord = MaxLocalVertex_Donor*nDim;
  nBuffer_Point = MaxLocalVertex_Donor;

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  if (faces) {
    SU2_MPI::Allgather(Buffer_Send_Normal, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Normal, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  }
#else
  for (iVertex = 0; iVertex < nBuffer_Coord; iVertex++)
    Buffer_Receive_Coord[iVertex] = Buffer_Send_Coord[iVertex];

  for (iVertex = 0; iVertex < nBuffer_Point; iVertex++)
    Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];

  if (faces) {
    for (iVertex = 0; iVertex < nBuffer_Coord; iVertex++)
      Buffer_Receive_Normal[iVertex] = Buffer_Send_Normal[iVertex];
  }
#endif
}

int CInterpolator::Find_InterfaceMarker(CConfig *config, unsigned short val_marker_interface) {
    
  unsigned short nMarker = config->GetnMarker_All();
  unsigned short iMarker;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    /*--- If the tag GetMarker_All_ZoneInterface(iMarker) equals the index we are looping at ---*/
    if (config->GetMarker_All_ZoneInterface(iMarker) == val_marker_interface ) {

      /*--- We have identified the identifier for the interface marker ---*/
      return iMarker;
    }
  }

  return -1;
}


void CInterpolator::ReconstructBoundary(unsigned long val_zone, int val_marker){
    
  CGeometry *geom = Geometry[val_zone][MESH_0];
    
  unsigned long iVertex, jVertex, kVertex;
    
  unsigned long count, iTmp, *uptr, dPoint, EdgeIndex, jEdge, nEdges, nNodes, nVertex, iDim, nDim, iPoint;
   
  unsigned long nGlobalLinkedNodes, nLocalVertex, nLocalLinkedNodes;
  
  nDim = geom->GetnDim();
  
  if( val_marker != -1 )
    nVertex  = geom->GetnVertex(  val_marker  );
  else
    nVertex  = 0;
      
    
  su2double *Buffer_Send_Coord           = new su2double     [ nVertex * nDim ];
  unsigned long *Buffer_Send_GlobalPoint = new unsigned long [ nVertex ];
  
  unsigned long *Buffer_Send_nLinkedNodes       = new unsigned long [ nVertex ];
  unsigned long *Buffer_Send_StartLinkedNodes   = new unsigned long [ nVertex ];
  unsigned long **Aux_Send_Map                  = new unsigned long*[ nVertex ];

#ifdef HAVE_MPI
  int nProcessor = size, iRank;
  unsigned long iTmp2, tmp_index, tmp_index_2;
#endif
        
  /*--- Copy coordinates and point to the auxiliar vector ---*/
  
  nGlobalVertex     = 0;
  nLocalVertex      = 0;
  nLocalLinkedNodes = 0;
  
  for (iVertex = 0; iVertex < nVertex; iVertex++) {
    
    Buffer_Send_nLinkedNodes[iVertex] = 0;
    Aux_Send_Map[iVertex]             = NULL;
    
    iPoint = geom->vertex[val_marker][iVertex]->GetNode();
    
    if (geom->node[iPoint]->GetDomain()) {
      Buffer_Send_GlobalPoint[nLocalVertex] = geom->node[iPoint]->GetGlobalIndex();
      
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[nLocalVertex*nDim+iDim] = geom->node[iPoint]->GetCoord(iDim);
   
      nNodes = 0;
      nEdges = geom->node[iPoint]->GetnPoint();
        
      for (jEdge = 0; jEdge < nEdges; jEdge++){
        EdgeIndex = geom->node[iPoint]->GetEdge(jEdge);

        if( iPoint == geom->edge[EdgeIndex]->GetNode(0) )
          dPoint = geom->edge[EdgeIndex]->GetNode(1);
        else
          dPoint = geom->edge[EdgeIndex]->GetNode(0);

        if ( geom->node[dPoint]->GetVertex(val_marker) != -1 )
          nNodes++;
      }

      Buffer_Send_StartLinkedNodes[nLocalVertex] = nLocalLinkedNodes;
      Buffer_Send_nLinkedNodes[nLocalVertex]     = nNodes;

      nLocalLinkedNodes += nNodes;

      Aux_Send_Map[nLocalVertex] = new unsigned long[ nNodes ];
      nNodes = 0;

      for (jEdge = 0; jEdge < nEdges; jEdge++){    
        EdgeIndex = geom->node[iPoint]->GetEdge(jEdge);

        if( iPoint == geom->edge[EdgeIndex]->GetNode(0) )
          dPoint = geom->edge[EdgeIndex]->GetNode(1);
        else
          dPoint = geom->edge[EdgeIndex]->GetNode(0);                

        if ( geom->node[dPoint]->GetVertex(val_marker) != -1 ){    
          Aux_Send_Map[nLocalVertex][nNodes] = geom->node[dPoint]->GetGlobalIndex();
          nNodes++;
        }
      }  
      nLocalVertex++;
    }
  }
    
  unsigned long *Buffer_Send_LinkedNodes = new unsigned long [ nLocalLinkedNodes ];

  nLocalLinkedNodes = 0;

  for (iVertex = 0; iVertex < nLocalVertex; iVertex++){
    for (jEdge = 0; jEdge < Buffer_Send_nLinkedNodes[iVertex]; jEdge++){
      Buffer_Send_LinkedNodes[nLocalLinkedNodes] = Aux_Send_Map[iVertex][jEdge];
      nLocalLinkedNodes++;
    }
  }
    
 for (iVertex = 0; iVertex < nVertex; iVertex++){
    if( Aux_Send_Map[iVertex] != NULL )
      delete [] Aux_Send_Map[iVertex];
  }
  delete [] Aux_Send_Map; Aux_Send_Map = NULL;

  /*--- Reconstruct  boundary by gathering data from all ranks ---*/

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(     &nLocalVertex,      &nGlobalVertex, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalLinkedNodes, &nGlobalLinkedNodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nGlobalVertex      = nLocalVertex;
  nGlobalLinkedNodes = nLocalLinkedNodes;
#endif 

  Buffer_Receive_Coord       = new su2double    [ nGlobalVertex * nDim ];
  Buffer_Receive_GlobalPoint = new unsigned long[ nGlobalVertex ];
  Buffer_Receive_Proc        = new unsigned long[ nGlobalVertex ];
   
  Buffer_Receive_nLinkedNodes     = new unsigned long[ nGlobalVertex ];
  Buffer_Receive_LinkedNodes      = new unsigned long[ nGlobalLinkedNodes   ];
  Buffer_Receive_StartLinkedNodes = new unsigned long[ nGlobalVertex ];

#ifdef HAVE_MPI
  if (rank == MASTER_NODE){

    for (iVertex = 0; iVertex < nDim*nLocalVertex; iVertex++)
      Buffer_Receive_Coord[iVertex]  = Buffer_Send_Coord[iVertex];

    for (iVertex = 0; iVertex < nLocalVertex; iVertex++){
      Buffer_Receive_GlobalPoint[iVertex]      = Buffer_Send_GlobalPoint[iVertex];
      Buffer_Receive_Proc[iVertex]             = MASTER_NODE;
      Buffer_Receive_nLinkedNodes[iVertex]     = Buffer_Send_nLinkedNodes[iVertex];
      Buffer_Receive_StartLinkedNodes[iVertex] = Buffer_Send_StartLinkedNodes[iVertex];
    }
      
    for (iVertex = 0; iVertex < nLocalLinkedNodes; iVertex++)
      Buffer_Receive_LinkedNodes[iVertex] = Buffer_Send_LinkedNodes[iVertex];
 
    tmp_index   = nLocalVertex;
    tmp_index_2 = nLocalLinkedNodes;

    for(iRank = 1; iRank < nProcessor; iRank++){
       
      SU2_MPI::Recv(                           &iTmp2,     1, MPI_UNSIGNED_LONG, iRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_LinkedNodes[tmp_index_2], iTmp2, MPI_UNSIGNED_LONG, iRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      SU2_MPI::Recv(                         &iTmp,         1, MPI_UNSIGNED_LONG, iRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_Coord[tmp_index*nDim], nDim*iTmp,        MPI_DOUBLE, iRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      SU2_MPI::Recv(     &Buffer_Receive_GlobalPoint[tmp_index], iTmp, MPI_UNSIGNED_LONG, iRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      SU2_MPI::Recv(    &Buffer_Receive_nLinkedNodes[tmp_index], iTmp, MPI_UNSIGNED_LONG, iRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      SU2_MPI::Recv(&Buffer_Receive_StartLinkedNodes[tmp_index], iTmp, MPI_UNSIGNED_LONG, iRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      for (iVertex = 0; iVertex < iTmp; iVertex++){
        Buffer_Receive_Proc[ tmp_index + iVertex ] = iRank;
        Buffer_Receive_StartLinkedNodes[ tmp_index + iVertex ] += tmp_index_2;
      }
        
      tmp_index   += iTmp;
      tmp_index_2 += iTmp2;
    }
  }
  else{
    SU2_MPI::Send(     &nLocalLinkedNodes,                 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
    SU2_MPI::Send(Buffer_Send_LinkedNodes, nLocalLinkedNodes, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
    
    SU2_MPI::Send(    &nLocalVertex,                   1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
    SU2_MPI::Send(Buffer_Send_Coord, nDim * nLocalVertex,        MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      
    SU2_MPI::Send(     Buffer_Send_GlobalPoint, nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
    SU2_MPI::Send(    Buffer_Send_nLinkedNodes, nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
    SU2_MPI::Send(Buffer_Send_StartLinkedNodes, nLocalVertex, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
  }    
#else
  for (iVertex = 0; iVertex < nDim * nGlobalVertex; iVertex++)
    Buffer_Receive_Coord[iVertex] = Buffer_Send_Coord[iVertex];
     
  for (iVertex = 0; iVertex < nGlobalVertex; iVertex++){
    Buffer_Receive_GlobalPoint[iVertex]      = Buffer_Send_GlobalPoint[iVertex];
    Buffer_Receive_Proc[iVertex]             = MASTER_NODE;
    Buffer_Receive_nLinkedNodes[iVertex]     = Buffer_Send_nLinkedNodes[iVertex];
    Buffer_Receive_StartLinkedNodes[iVertex] = Buffer_Send_StartLinkedNodes[iVertex];
  }
    
  for (iVertex = 0; iVertex < nGlobalLinkedNodes; iVertex++)
    Buffer_Receive_LinkedNodes[iVertex] = Buffer_Send_LinkedNodes[iVertex];
#endif 

  if (rank == MASTER_NODE){
    for (iVertex = 0; iVertex < nGlobalVertex; iVertex++){
      count = 0;
      uptr = &Buffer_Receive_LinkedNodes[ Buffer_Receive_StartLinkedNodes[iVertex] ];
      
      for (jVertex = 0; jVertex < Buffer_Receive_nLinkedNodes[iVertex]; jVertex++){
        iTmp = uptr[ jVertex ];
        for (kVertex = 0; kVertex < nGlobalVertex; kVertex++){
          if( Buffer_Receive_GlobalPoint[kVertex] == iTmp ){
            uptr[ jVertex ] = kVertex;
            count++;
            break;
          }
        }
          
        if( count != (jVertex+1) ){
          for (kVertex = jVertex; kVertex < Buffer_Receive_nLinkedNodes[iVertex]-1; kVertex++){
            uptr[ kVertex ] = uptr[ kVertex + 1];
          }
          Buffer_Receive_nLinkedNodes[iVertex]--;
          jVertex--;   
        }
      }
    }
  }

#ifdef HAVE_MPI    
  SU2_MPI::Bcast(      Buffer_Receive_Coord, nGlobalVertex * nDim,        MPI_DOUBLE, 0, MPI_COMM_WORLD);
  SU2_MPI::Bcast(Buffer_Receive_GlobalPoint, nGlobalVertex,        MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  SU2_MPI::Bcast(      Buffer_Receive_Proc, nGlobalVertex,        MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
  
  SU2_MPI::Bcast(    Buffer_Receive_nLinkedNodes,      nGlobalVertex, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  SU2_MPI::Bcast(Buffer_Receive_StartLinkedNodes,      nGlobalVertex, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  SU2_MPI::Bcast(     Buffer_Receive_LinkedNodes, nGlobalLinkedNodes, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
#endif
  
  if( Buffer_Send_Coord              != NULL) {delete [] Buffer_Send_Coord;            Buffer_Send_Coord            = NULL;} 
  if( Buffer_Send_GlobalPoint        != NULL) {delete [] Buffer_Send_GlobalPoint;      Buffer_Send_GlobalPoint      = NULL;}
  if( Buffer_Send_LinkedNodes        != NULL) {delete [] Buffer_Send_LinkedNodes;      Buffer_Send_LinkedNodes      = NULL;}
  if( Buffer_Send_nLinkedNodes       != NULL) {delete [] Buffer_Send_nLinkedNodes;     Buffer_Send_nLinkedNodes     = NULL;}
  if( Buffer_Send_StartLinkedNodes   != NULL) {delete [] Buffer_Send_StartLinkedNodes; Buffer_Send_StartLinkedNodes = NULL;}
}

bool CInterpolator::CheckInterfaceBoundary(int markDonor, int markTarget){
  
  int Donor_check, Target_check;
  
  #ifdef HAVE_MPI
    
  int *Buffer_Recv_mark = NULL;
  int iRank, nProcessor = size;
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[nProcessor];

  Donor_check  = -1;
  Target_check = -1;

  /*--- We gather a vector in MASTER_NODE to determine whether the boundary is not on the processor because of the partition or because the zone does not include it ---*/

  SU2_MPI::Gather(&markDonor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

  if (rank == MASTER_NODE)
    for (iRank = 0; iRank < nProcessor; iRank++)
      if( Buffer_Recv_mark[iRank] != -1 ){
        Donor_check = Buffer_Recv_mark[iRank];
        break;
      }

  SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);


  SU2_MPI::Gather(&markTarget, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

  if (rank == MASTER_NODE)
    for (iRank = 0; iRank < nProcessor; iRank++)
      if( Buffer_Recv_mark[iRank] != -1 ){
        Target_check = Buffer_Recv_mark[iRank];
        break;
      }


  SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
  
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;

#else
  Donor_check  = markDonor;
  Target_check = markTarget;
#endif

  if(Target_check == -1 || Donor_check == -1)
    return false;
  else 
    return true;
}

su2double CInterpolator::PointsDistance(su2double *point_i, su2double *point_j){

  /*--- Compute distance between 2 points ---*/

  unsigned short iDim, nDim = donor_geometry->GetnDim();
  su2double m;

  m = 0 ;
  for(iDim = 0; iDim < nDim; iDim++)
    m += (point_j[iDim] - point_i[iDim])*(point_j[iDim] - point_i[iDim]);

  return sqrt(m);
}

/* Nearest Neighbor Interpolator */
CNearestNeighbor::CNearestNeighbor(void):  CInterpolator() { }

CNearestNeighbor::CNearestNeighbor(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone) {

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);
}

CNearestNeighbor::~CNearestNeighbor() {}

void CNearestNeighbor::Set_TransferCoeff(CConfig **config) {

  int iProcessor, pProcessor, nProcessor = size;
  int markDonor, markTarget;

  unsigned short nDim, iMarkerInt, nMarkerInt, iDonor;    

  unsigned long nVertexDonor, nVertexTarget, Point_Target, jVertex, iVertexTarget;
  unsigned long Global_Point_Donor, pGlobalPoint=0;

  su2double *Coord_i, *Coord_j, dist, mindist, maxdist;

  /*--- Initialize variables --- */
  
  nMarkerInt = (int) ( config[donorZone]->GetMarker_n_ZoneInterface() / 2 );
  
  nDim = donor_geometry->GetnDim();

  iDonor = 0;
  
  Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];


  /*--- Cycle over nMarkersInt interface to determine communication pattern ---*/

  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {


    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);
      
    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if( !CheckInterfaceBoundary(markDonor, markTarget) )
      continue;

    if(markDonor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( markDonor );
    else
      nVertexDonor  = 0;
    
    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;
    
    Buffer_Send_nVertex_Donor  = new unsigned long [ 1 ];

    /* Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor */
    Determine_ArraySize(false, markDonor, markTarget, nVertexDonor, nDim);

    Buffer_Send_Coord          = new su2double     [ MaxLocalVertex_Donor * nDim ];
    Buffer_Send_GlobalPoint    = new unsigned long [ MaxLocalVertex_Donor ];
    Buffer_Receive_Coord       = new su2double     [ nProcessor * MaxLocalVertex_Donor * nDim ];
    Buffer_Receive_GlobalPoint = new unsigned long [ nProcessor * MaxLocalVertex_Donor ];

    /*-- Collect coordinates, global points, and normal vectors ---*/
    Collect_VertexInfo( false, markDonor, markTarget, nVertexDonor, nDim );

    /*--- Compute the closest point to a Near-Field boundary point ---*/
    maxdist = 0.0;

    for (iVertexTarget = 0; iVertexTarget < nVertexTarget; iVertexTarget++) {

      Point_Target = target_geometry->vertex[markTarget][iVertexTarget]->GetNode();

      if ( target_geometry->node[Point_Target]->GetDomain() ) {

        target_geometry->vertex[markTarget][iVertexTarget]->SetnDonorPoints(1);
        target_geometry->vertex[markTarget][iVertexTarget]->Allocate_DonorInfo(); // Possible meme leak?

        /*--- Coordinates of the boundary point ---*/
        Coord_i = target_geometry->node[Point_Target]->GetCoord();

        mindist    = 1E6; 
        pProcessor = 0;

        /*--- Loop over all the boundaries to find the pair ---*/

        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (jVertex = 0; jVertex < MaxLocalVertex_Donor; jVertex++) {
            Global_Point_Donor = iProcessor*MaxLocalVertex_Donor+jVertex;

            Coord_j = &Buffer_Receive_Coord[ Global_Point_Donor*nDim];

            dist = PointsDistance(Coord_i, Coord_j);

            if (dist < mindist) {
              mindist = dist; pProcessor = iProcessor; pGlobalPoint = Buffer_Receive_GlobalPoint[Global_Point_Donor];
            }

            if (dist == 0.0) break;
          }

        }

        /*--- Store the value of the pair ---*/
        maxdist = max(maxdist, mindist);
        target_geometry->vertex[markTarget][iVertexTarget]->SetInterpDonorPoint(iDonor, pGlobalPoint);
        target_geometry->vertex[markTarget][iVertexTarget]->SetInterpDonorProcessor(iDonor, pProcessor);
        target_geometry->vertex[markTarget][iVertexTarget]->SetDonorCoeff(iDonor, 1.0);
      }
    }

    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_GlobalPoint;
    
    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_GlobalPoint;

    delete[] Buffer_Send_nVertex_Donor;

  }

  delete[] Buffer_Receive_nVertex_Donor;
}



CIsoparametric::CIsoparametric(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone)  :  CInterpolator(geometry_container, config, iZone, jZone) {

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

  /*--- For fluid-structure interaction data interpolated with have nDim dimensions ---*/
 // InitializeData(Zones,nDim);
}

CIsoparametric::~CIsoparametric() {}

void CIsoparametric::Set_TransferCoeff(CConfig **config) {
  unsigned long iVertex, jVertex;
  unsigned long  dPoint, inode, jElem, nElem;
  unsigned short iDim, iDonor=0, iFace;

  unsigned short nDim = donor_geometry->GetnDim();

  unsigned short nMarkerInt;
  unsigned short iMarkerInt;

  int markDonor=0, markTarget=0;

  long donor_elem=0, temp_donor=0;
  unsigned int nNodes=0;
  /*--- Restricted to 2-zone for now ---*/
  unsigned int nFaces=1; //For 2D cases, we want to look at edges, not faces, as the 'interface'
  bool face_on_marker=true;

  unsigned long nVertexDonor = 0, nVertexTarget= 0;
  unsigned long Point_Target = 0;

  unsigned long iVertexDonor, iPointDonor = 0;
  unsigned long jGlobalPoint = 0;
  int iProcessor;

  unsigned long nLocalFace_Donor = 0, nLocalFaceNodes_Donor=0;

  unsigned long faceindex;

  su2double dist = 0.0, mindist=1E6, *Coord, *Coord_i;
  su2double myCoeff[10]; // Maximum # of donor points
  su2double  *Normal;
  su2double *projected_point = new su2double[nDim];
  su2double tmp, tmp2;
  su2double storeCoeff[10];
  unsigned long storeGlobal[10];
  int storeProc[10];

  int nProcessor = size;
  Coord = new su2double[nDim];
  Normal = new su2double[nDim];

  nMarkerInt = (config[donorZone]->GetMarker_n_ZoneInterface())/2;

  /*--- For the number of markers on the interface... ---*/
  for (iMarkerInt=1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
    /*--- Procedure:
    * -Loop through vertices of the aero grid
    * -Find nearest element and allocate enough space in the aero grid donor point info
    *    -set the transfer coefficient values
    */

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);
      
    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if( !CheckInterfaceBoundary(markDonor, markTarget) )
      continue;

    if(markDonor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( markDonor );
    else
      nVertexDonor  = 0;

    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;
    
    Buffer_Send_nVertex_Donor    = new unsigned long [1];
    Buffer_Send_nFace_Donor      = new unsigned long [1];
    Buffer_Send_nFaceNodes_Donor = new unsigned long [1];

    Buffer_Receive_nVertex_Donor    = new unsigned long [nProcessor];
    Buffer_Receive_nFace_Donor      = new unsigned long [nProcessor];
    Buffer_Receive_nFaceNodes_Donor = new unsigned long [nProcessor];

    /* Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor */
    Determine_ArraySize(true, markDonor, markTarget, nVertexDonor, nDim);

    Buffer_Send_Coord       = new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_Normal      = new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_GlobalPoint = new unsigned long [MaxLocalVertex_Donor];

    Buffer_Receive_Coord       = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_Normal      = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Donor];

    /*-- Collect coordinates, global points, and normal vectors ---*/
    Collect_VertexInfo(true, markDonor,markTarget,nVertexDonor,nDim);

    Buffer_Send_FaceIndex    = new unsigned long[MaxFace_Donor];
    Buffer_Send_FaceNodes    = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_FaceProc     = new unsigned long[MaxFaceNodes_Donor];

    Buffer_Receive_FaceIndex = new unsigned long[MaxFace_Donor*nProcessor];
    Buffer_Receive_FaceNodes = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_FaceProc  = new unsigned long[MaxFaceNodes_Donor*nProcessor];

    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;

    /*--- Collect Face info ---*/

    for (iVertex = 0; iVertex < MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex] = 0;
    }
    for (iVertex=0; iVertex<MaxFaceNodes_Donor; iVertex++) {
      Buffer_Send_FaceNodes[iVertex] = 0;
      Buffer_Send_FaceProc[iVertex]  = 0;
    }

    Buffer_Send_FaceIndex[0] = rank * MaxFaceNodes_Donor;

    if (nDim==2) nNodes=2;

    for (iVertexDonor = 0; iVertexDonor < nVertexDonor; iVertexDonor++) {
      iPointDonor = donor_geometry->vertex[markDonor][iVertexDonor]->GetNode();

      if (donor_geometry->node[iPointDonor]->GetDomain()) {

    if (nDim==3)  nElem = donor_geometry->node[iPointDonor]->GetnElem();
    else          nElem =donor_geometry->node[iPointDonor]->GetnPoint();

    for (jElem=0; jElem < nElem; jElem++) {
      if (nDim==3) {
        temp_donor = donor_geometry->node[iPointDonor]->GetElem(jElem);
        nFaces = donor_geometry->elem[temp_donor]->GetnFaces();
        for (iFace=0; iFace<nFaces; iFace++) {
          /*-- Determine whether this face/edge is on the marker --*/
          face_on_marker=true;
          nNodes = donor_geometry->elem[temp_donor]->GetnNodesFace(iFace);
          for (iDonor=0; iDonor<nNodes; iDonor++) {
            inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, iDonor);
            dPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
            face_on_marker = (face_on_marker && (donor_geometry->node[dPoint]->GetVertex(markDonor) !=-1));
          }

          if (face_on_marker ) {
            for (iDonor=0; iDonor<nNodes; iDonor++) {
              inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, iDonor);
              dPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
              // Match node on the face to the correct global index
              jGlobalPoint=donor_geometry->node[dPoint]->GetGlobalIndex();
              for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                for (jVertex = 0; jVertex < Buffer_Receive_nVertex_Donor[iProcessor]; jVertex++) {
                  if (jGlobalPoint ==Buffer_Receive_GlobalPoint[MaxLocalVertex_Donor*iProcessor+jVertex]) {
                    Buffer_Send_FaceNodes[nLocalFaceNodes_Donor]=MaxLocalVertex_Donor*iProcessor+jVertex;
                    Buffer_Send_FaceProc[nLocalFaceNodes_Donor]=iProcessor;
                  }
                }
              }
              nLocalFaceNodes_Donor++; // Increment total number of face-nodes / processor
            }
            /* Store the indices */
            Buffer_Send_FaceIndex[nLocalFace_Donor+1] = Buffer_Send_FaceIndex[nLocalFace_Donor]+nNodes;
            nLocalFace_Donor++; // Increment number of faces / processor
          }
        }
      }
      else {
        /*-- Determine whether this face/edge is on the marker --*/
        face_on_marker=true;
        for (iDonor=0; iDonor<nNodes; iDonor++) {
          inode = donor_geometry->node[iPointDonor]->GetEdge(jElem);
          dPoint = donor_geometry->edge[inode]->GetNode(iDonor);
          face_on_marker = (face_on_marker && (donor_geometry->node[dPoint]->GetVertex(markDonor) !=-1));
        }
        if (face_on_marker ) {
          for (iDonor=0; iDonor<nNodes; iDonor++) {
            inode = donor_geometry->node[iPointDonor]->GetEdge(jElem);
            dPoint = donor_geometry->edge[inode]->GetNode(iDonor);
            // Match node on the face to the correct global index
            jGlobalPoint=donor_geometry->node[dPoint]->GetGlobalIndex();
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
              for (jVertex = 0; jVertex < Buffer_Receive_nVertex_Donor[iProcessor]; jVertex++) {
                if (jGlobalPoint ==Buffer_Receive_GlobalPoint[MaxLocalVertex_Donor*iProcessor+jVertex]) {
                  Buffer_Send_FaceNodes[nLocalFaceNodes_Donor]=MaxLocalVertex_Donor*iProcessor+jVertex;
                  Buffer_Send_FaceProc[nLocalFaceNodes_Donor]=iProcessor;
                }
              }
            }
            nLocalFaceNodes_Donor++; // Increment total number of face-nodes / processor
          }
          /* Store the indices */
          Buffer_Send_FaceIndex[nLocalFace_Donor+1] = Buffer_Send_FaceIndex[nLocalFace_Donor]+nNodes;
          nLocalFace_Donor++; // Increment number of faces / processor
        }
      }
    }
      }
    }

    //Buffer_Send_FaceIndex[nLocalFace_Donor+1] = MaxFaceNodes_Donor*rank+nLocalFaceNodes_Donor;
#ifdef HAVE_MPI
    SU2_MPI::Allgather(Buffer_Send_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceProc, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceProc, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
    for (iFace=0; iFace<MaxFace_Donor; iFace++) {
      Buffer_Receive_FaceIndex[iFace] = Buffer_Send_FaceIndex[iFace];
    }
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++)
      Buffer_Receive_FaceNodes[iVertex] = Buffer_Send_FaceNodes[iVertex];
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++)
      Buffer_Receive_FaceProc[iVertex] = Buffer_Send_FaceProc[iVertex];
#endif

    /*--- Loop over the vertices on the target Marker ---*/
    for (iVertex = 0; iVertex<nVertexTarget; iVertex++) {
      mindist=1E6;
      for (unsigned short iCoeff=0; iCoeff<10; iCoeff++) {
    storeCoeff[iCoeff]=0;
      }
      Point_Target = target_geometry->vertex[markTarget][iVertex]->GetNode();

      if (target_geometry->node[Point_Target]->GetDomain()) {

    Coord_i = target_geometry->node[Point_Target]->GetCoord();
    /*---Loop over the faces previously communicated/stored ---*/
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {

      nFaces = (unsigned int)Buffer_Receive_nFace_Donor[iProcessor];

      for (iFace = 0; iFace< nFaces; iFace++) {
        /*--- ---*/

        nNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1] -
                (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace];

        su2double *X = new su2double[nNodes*nDim];
        faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
        for (iDonor=0; iDonor<nNodes; iDonor++) {
          jVertex = Buffer_Receive_FaceNodes[iDonor+faceindex]; // index which points to the stored coordinates, global points
          for (iDim=0; iDim<nDim; iDim++) {
            X[iDim*nNodes+iDonor]=
                Buffer_Receive_Coord[jVertex*nDim+iDim];
          }
        }
        jVertex = Buffer_Receive_FaceNodes[faceindex];

        for (iDim=0; iDim<nDim; iDim++) {
          Normal[iDim] = Buffer_Receive_Normal[jVertex*nDim+iDim];
        }

        /* Project point used for case where surfaces are not exactly coincident, where
         * the point is assumed connected by a rigid rod normal to the surface.
         */
        tmp = 0;
        tmp2=0;
        for (iDim=0; iDim<nDim; iDim++) {
          tmp+=Normal[iDim]*Normal[iDim];
          tmp2+=Normal[iDim]*(Coord_i[iDim]-X[iDim*nNodes]);
        }
        tmp = 1/tmp;
        tmp2 = tmp2*sqrt(tmp);
        for (iDim=0; iDim<nDim; iDim++) {
          // projection of \vec{q} onto plane defined by \vec{n} and \vec{p}:
          // \vec{q} - \vec{n} ( (\vec{q}-\vec{p} ) \cdot \vec{n})
          // tmp2 = ( (\vec{q}-\vec{p} ) \cdot \vec{N})
          // \vec{n} = \vec{N}/(|N|), tmp = 1/|N|^2
          projected_point[iDim]=Coord_i[iDim] + Normal[iDim]*tmp2*tmp;
        }

        Isoparameters(nDim, nNodes, X, projected_point,myCoeff);

        /*--- Find distance to the interpolated point ---*/
        dist = 0.0;
        for (iDim=0; iDim<nDim; iDim++) {
          Coord[iDim] = Coord_i[iDim];
          for(iDonor=0; iDonor< nNodes; iDonor++) {
            Coord[iDim]-=myCoeff[iDonor]*X[iDim*nNodes+iDonor];
          }
          dist+=pow(Coord[iDim],2.0);
        }

        /*--- If the dist is shorter than last closest (and nonzero nodes are on the boundary), update ---*/
        if (dist<mindist ) {
          /*--- update last dist ---*/
          mindist = dist;
          /*--- Store info ---*/
          donor_elem = temp_donor;
          target_geometry->vertex[markTarget][iVertex]->SetDonorElem(donor_elem); // in 2D is nearest neighbor
          target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);
          for (iDonor=0; iDonor<nNodes; iDonor++) {
            storeCoeff[iDonor] = myCoeff[iDonor];
            jVertex = Buffer_Receive_FaceNodes[faceindex+iDonor];
            storeGlobal[iDonor] =Buffer_Receive_GlobalPoint[jVertex];
            storeProc[iDonor] = (int)Buffer_Receive_FaceProc[faceindex+iDonor];
          }
        }
      
        delete [] X;
      }
    }
    /*--- Set the appropriate amount of memory and fill ---*/
    nNodes =target_geometry->vertex[markTarget][iVertex]->GetnDonorPoints();
    target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

    for (iDonor=0; iDonor<nNodes; iDonor++) {
      target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint(iDonor,storeGlobal[iDonor]);
      //cout <<rank << " Global Point " << Global_Point<<" iDonor " << iDonor <<" coeff " << coeff <<" gp " << pGlobalPoint << endl;
      target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(iDonor,storeCoeff[iDonor]);
      target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, storeProc[iDonor]);
    }
      }
    }

    delete[] Buffer_Send_nVertex_Donor;
    delete[] Buffer_Send_nFace_Donor;
    delete[] Buffer_Send_nFaceNodes_Donor;

    delete[] Buffer_Receive_nVertex_Donor;
    delete[] Buffer_Receive_nFace_Donor;
    delete[] Buffer_Receive_nFaceNodes_Donor;

    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_Normal;
    delete[] Buffer_Send_GlobalPoint;

    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_Normal;
    delete[] Buffer_Receive_GlobalPoint;

    delete[] Buffer_Send_FaceIndex;
    delete[] Buffer_Send_FaceNodes;
    delete[] Buffer_Send_FaceProc;

    delete[] Buffer_Receive_FaceIndex;
    delete[] Buffer_Receive_FaceNodes;
    delete[] Buffer_Receive_FaceProc;
  }
  delete [] Coord;
  delete [] Normal;
  
  delete [] projected_point;
}

void CIsoparametric::Isoparameters(unsigned short nDim, unsigned short nDonor,
    su2double *X, su2double *xj, su2double *isoparams) {
  short iDonor,iDim,k; // indices
  su2double tmp, tmp2;
  
  su2double *x     = new su2double[nDim+1];
  su2double *x_tmp = new su2double[nDim+1];
  su2double *Q     = new su2double[nDonor*nDonor];
  su2double *R     = new su2double[nDonor*nDonor];
  su2double *A     = new su2double[nDim+1*nDonor];
  su2double *A2    = NULL;
  su2double *x2    = new su2double[nDim+1];
  
  bool *test  = new bool[nDim+1];
  bool *testi = new bool[nDim+1];
  
  su2double eps = 1E-10;
  
  short n = nDim+1;

  if (nDonor>2) {
    /*--- Create Matrix A: 1st row all 1's, 2nd row x coordinates, 3rd row y coordinates, etc ---*/
    /*--- Right hand side is [1, \vec{x}']'---*/
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      isoparams[iDonor]=0;
      A[iDonor] = 1.0;
      for (iDim=0; iDim<n; iDim++)
        A[(iDim+1)*nDonor+iDonor]=X[iDim*nDonor+iDonor];
    }

    x[0] = 1.0;
    for (iDim=0; iDim<nDim; iDim++)
      x[iDim+1]=xj[iDim];

    /*--- Eliminate degenerate rows:
     * for example, if z constant including the z values will make the system degenerate
     * TODO: improve efficiency of this loop---*/
    test[0]=true; // always keep the 1st row
    for (iDim=1; iDim<nDim+1; iDim++) {
      // Test this row against all previous
      test[iDim]=true; // Assume that it is not degenerate
      for (k=0; k<iDim; k++) {
        tmp=0; tmp2=0;
        for (iDonor=0;iDonor<nDonor;iDonor++) {
          tmp+= A[iDim*nDonor+iDonor]*A[iDim*nDonor+iDonor];
          tmp2+=A[k*nDonor+iDonor]*A[k*nDonor+iDonor];
        }
        tmp  = pow(tmp,0.5);
        tmp2 = pow(tmp2,0.5);
        testi[k]=false;
        for (iDonor=0; iDonor<nDonor; iDonor++) {
          // If at least one ratio is non-matching row iDim is not degenerate w/ row k
          if (A[iDim*nDonor+iDonor]/tmp != A[k*nDonor+iDonor]/tmp2)
            testi[k]=true;
        }
        // If any of testi (k<iDim) are false, row iDim is degenerate
        test[iDim]=(test[iDim] && testi[k]);
      }
      if (!test[iDim]) n--;
    }

    /*--- Initialize A2 now that we might have a smaller system --*/
    A2 = new su2double[n*nDonor];
    iDim=0;
    /*--- Copy only the rows that are non-degenerate ---*/
    for (k=0; k<nDim+1; k++) {
      if (test[k]) {
        for (iDonor=0;iDonor<nDonor;iDonor++ ) {
          A2[nDonor*iDim+iDonor]=A[nDonor*k+iDonor];
        }
        x2[iDim]=x[k];
        iDim++;
      }
    }
    /*--- Initialize Q,R to 0 --*/
    for (k=0; k<nDonor*nDonor; k++) {
      Q[k]=0;
      R[k]=0;
    }
    /*--- TODO: make this loop more efficient ---*/
    /*--- Solve for rectangular Q1 R1 ---*/
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      tmp=0;
      for (iDim=0; iDim<n; iDim++)
        tmp += (A2[iDim*nDonor+iDonor])*(A2[iDim*nDonor+iDonor]);

      R[iDonor*nDonor+iDonor]= pow(tmp,0.5);
      if (tmp>eps && iDonor<n) {
        for (iDim=0; iDim<n; iDim++)
          Q[iDim*nDonor+iDonor]=A2[iDim*nDonor+iDonor]/R[iDonor*nDonor+iDonor];
      }
      else if (tmp!=0) {
        for (iDim=0; iDim<n; iDim++)
          Q[iDim*nDonor+iDonor]=A2[iDim*nDonor+iDonor]/tmp;
      }
      for (iDim=iDonor+1; iDim<nDonor; iDim++) {
        tmp=0;
        for (k=0; k<n; k++)
          tmp+=A2[k*nDonor+iDim]*Q[k*nDonor+iDonor];

        R[iDonor*nDonor+iDim]=tmp;

        for (k=0; k<n; k++)
          A2[k*nDonor+iDim]=A2[k*nDonor+iDim]-Q[k*nDonor+iDonor]*R[iDonor*nDonor+iDim];
      }
    }
    /*--- x_tmp = Q^T * x2 ---*/
    for (iDonor=0; iDonor<nDonor; iDonor++)
      x_tmp[iDonor]=0.0;
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      for (iDim=0; iDim<n; iDim++)
        x_tmp[iDonor]+=Q[iDim*nDonor+iDonor]*x2[iDim];
    }

    /*--- solve x_tmp = R*isoparams for isoparams: upper triangular system ---*/
    for (iDonor = n-1; iDonor>=0; iDonor--) {
      if (R[iDonor*nDonor+iDonor]>eps)
        isoparams[iDonor]=x_tmp[iDonor]/R[iDonor*nDonor+iDonor];
      else
        isoparams[iDonor]=0;
      for (k=0; k<iDonor; k++)
        x_tmp[k]=x_tmp[k]-R[k*nDonor+iDonor]*isoparams[iDonor];
    }
  }
  else {
    /*-- For 2-donors (lines) it is simpler: */
    tmp =  pow(X[0*nDonor+0]- X[0*nDonor+1],2.0);
    tmp += pow(X[1*nDonor+0]- X[1*nDonor+1],2.0);
    tmp = sqrt(tmp);

    tmp2 = pow(X[0*nDonor+0] - xj[0],2.0);
    tmp2 += pow(X[1*nDonor+0] - xj[1],2.0);
    tmp2 = sqrt(tmp2);
    isoparams[1] = tmp2/tmp;

    tmp2 = pow(X[0*nDonor+1] - xj[0],2.0);
    tmp2 += pow(X[1*nDonor+1] - xj[1],2.0);
    tmp2 = sqrt(tmp2);
    isoparams[0] = tmp2/tmp;
  }

  /*--- Isoparametric coefficients have been calculated. Run checks to eliminate outside-element issues ---*/
  if (nDonor==4) {
    /*-- Bilinear coordinates, bounded by [-1,1] ---*/
    su2double xi, eta;
    xi = (1.0-isoparams[0]/isoparams[1])/(1.0+isoparams[0]/isoparams[1]);
    eta = 1- isoparams[2]*4/(1+xi);
    if (xi>1.0) xi=1.0;
    if (xi<-1.0) xi=-1.0;
    if (eta>1.0) eta=1.0;
    if (eta<-1.0) eta=-1.0;
    isoparams[0]=0.25*(1-xi)*(1-eta);
    isoparams[1]=0.25*(1+xi)*(1-eta);
    isoparams[2]=0.25*(1+xi)*(1+eta);
    isoparams[3]=0.25*(1-xi)*(1+eta);

  }
  if (nDonor<4) {
    tmp = 0.0; // value for normalization
    tmp2=0; // check for maximum value, to be used to id nearest neighbor if necessary
    k=0; // index for maximum value
    for (iDonor=0; iDonor< nDonor; iDonor++) {
      if (isoparams[iDonor]>tmp2) {
        k=iDonor;
        tmp2=isoparams[iDonor];
      }
      // [0,1]
      if (isoparams[iDonor]<0) isoparams[iDonor]=0;
      if (isoparams[iDonor]>1) isoparams[iDonor] = 1;
      tmp +=isoparams[iDonor];
    }
    if (tmp>0)
      for (iDonor=0; iDonor< nDonor; iDonor++)
        isoparams[iDonor]=isoparams[iDonor]/tmp;
    else {
      isoparams[k] = 1.0;
    }
  }
  
  delete [] x;
  delete [] x_tmp;
  delete [] Q;
  delete [] R;
  delete [] A;
  if (A2 != NULL) delete [] A2;
  delete [] x2;
  
  delete [] test;
  delete [] testi;

}


/* Mirror Interpolator */
CMirror::CMirror(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone) {

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

}

CMirror::~CMirror() {}

void CMirror::Set_TransferCoeff(CConfig **config) {
  unsigned long iVertex, jVertex;
  unsigned long iPoint;
  unsigned short iDonor=0, iFace=0, iTarget=0;

  unsigned short nMarkerInt;
  unsigned short iMarkerInt;

  int markDonor=0, markTarget=0;

  unsigned int nNodes=0, iNodes=0;
  unsigned long nVertexDonor = 0, nVertexTarget= 0;
  unsigned long Point_Donor = 0;
  unsigned long Global_Point = 0;
  unsigned long pGlobalPoint = 0;
  int iProcessor;

  unsigned long nLocalFace_Donor = 0, nLocalFaceNodes_Donor=0;

  unsigned long faceindex;

  int nProcessor = size;

  su2double *Buffer_Send_Coeff, *Buffer_Receive_Coeff;
  su2double coeff;

  /*--- Number of markers on the interface ---*/
  nMarkerInt = (config[targetZone]->GetMarker_n_ZoneInterface())/2;

  /*--- For the number of markers on the interface... ---*/
  for (iMarkerInt=1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
   /*--- Procedure:
    * -Loop through vertices of the aero grid
    * -Find nearest element and allocate enough space in the aero grid donor point info
    *    -set the transfer coefficient values
    */

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if( !CheckInterfaceBoundary(markDonor, markTarget) )
      continue;

    if(markDonor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( markDonor );
    else
      nVertexDonor  = 0;

    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;

    /*-- Collect the number of donor nodes: re-use 'Face' containers --*/
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;
    for (jVertex = 0; jVertex<nVertexDonor; jVertex++) {
      Point_Donor =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex

      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        nLocalFaceNodes_Donor+=nNodes;
        nLocalFace_Donor++;
      }
    }
    Buffer_Send_nFace_Donor= new unsigned long [1];
    Buffer_Send_nFaceNodes_Donor= new unsigned long [1];

    Buffer_Receive_nFace_Donor = new unsigned long [nProcessor];
    Buffer_Receive_nFaceNodes_Donor = new unsigned long [nProcessor];

    Buffer_Send_nFace_Donor[0] = nLocalFace_Donor;
    Buffer_Send_nFaceNodes_Donor[0] = nLocalFaceNodes_Donor;

    /*--- Send Interface vertex information --*/
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &MaxFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor, &MaxFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFace_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MaxFace_Donor++;
#else
    nGlobalFace_Donor       = nLocalFace_Donor;
    nGlobalFaceNodes_Donor  = nLocalFaceNodes_Donor;
    MaxFaceNodes_Donor      = nLocalFaceNodes_Donor;
    MaxFace_Donor           = nLocalFace_Donor+1;
    Buffer_Receive_nFace_Donor[0] = Buffer_Send_nFace_Donor[0];
    Buffer_Receive_nFaceNodes_Donor[0] = Buffer_Send_nFaceNodes_Donor[0];
#endif

    /*-- Send donor info --*/
    Buffer_Send_FaceIndex   = new unsigned long[MaxFace_Donor];
    Buffer_Send_FaceNodes   = new unsigned long[MaxFaceNodes_Donor];
    //Buffer_Send_FaceProc    = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_GlobalPoint = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_Coeff       = new su2double[MaxFaceNodes_Donor];

    Buffer_Receive_FaceIndex= new unsigned long[MaxFace_Donor*nProcessor];
    Buffer_Receive_FaceNodes= new unsigned long[MaxFaceNodes_Donor*nProcessor];
    //Buffer_Receive_FaceProc = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_GlobalPoint = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_Coeff    = new su2double[MaxFaceNodes_Donor*nProcessor];

    for (iVertex=0; iVertex<MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex]=0;
    }
    for (iVertex=0; iVertex<MaxFaceNodes_Donor; iVertex++) {
      Buffer_Send_FaceNodes[iVertex]=0;
      //Buffer_Send_FaceProc[iVertex]=0;
      Buffer_Send_GlobalPoint[iVertex]=0;
      Buffer_Send_Coeff[iVertex]=0.0;
    }
    for (iVertex=0; iVertex<MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex]=0;
    }

    Buffer_Send_FaceIndex[0]=rank*MaxFaceNodes_Donor;
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;

    for (jVertex = 0; jVertex<nVertexDonor; jVertex++) {

      Point_Donor =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex
      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        for (iDonor=0; iDonor<nNodes; iDonor++) {
          Buffer_Send_FaceNodes[nLocalFaceNodes_Donor] = donor_geometry->node[Point_Donor]->GetGlobalIndex();
          Buffer_Send_GlobalPoint[nLocalFaceNodes_Donor] =
              donor_geometry->vertex[markDonor][jVertex]->GetInterpDonorPoint(iDonor);
          Buffer_Send_Coeff[nLocalFaceNodes_Donor] =
              donor_geometry->vertex[markDonor][jVertex]->GetDonorCoeff(iDonor);
          nLocalFaceNodes_Donor++;
        }
        Buffer_Send_FaceIndex[nLocalFace_Donor+1] =Buffer_Send_FaceIndex[nLocalFace_Donor]+nNodes;
        nLocalFace_Donor++;
      }
    }

#ifdef HAVE_MPI
    SU2_MPI::Allgather(Buffer_Send_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_GlobalPoint, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG,Buffer_Receive_GlobalPoint, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Coeff, MaxFaceNodes_Donor, MPI_DOUBLE,Buffer_Receive_Coeff, MaxFaceNodes_Donor, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
    for (iFace=0; iFace<MaxFace_Donor; iFace++) {
      Buffer_Receive_FaceIndex[iFace] = Buffer_Send_FaceIndex[iFace];
    }
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++) {
      Buffer_Receive_FaceNodes[iVertex] = Buffer_Send_FaceNodes[iVertex];
      Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];
      Buffer_Receive_Coeff[iVertex] = Buffer_Send_Coeff[iVertex];
    }
#endif
    /*--- Loop over the vertices on the target Marker ---*/
    for (iVertex = 0; iVertex<nVertexTarget; iVertex++) {

      iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();
      if (target_geometry->node[iPoint]->GetDomain()) {
        Global_Point = target_geometry->node[iPoint]->GetGlobalIndex();
        nNodes = 0;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iFace = 0; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {
            faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
            iNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1]- (unsigned int)faceindex;
            for (iTarget=0; iTarget<iNodes; iTarget++) {
              if (Global_Point == Buffer_Receive_GlobalPoint[faceindex+iTarget])
                nNodes++;
              //coeff =Buffer_Receive_Coeff[faceindex+iDonor];
            }
          }
        }

        target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);
        target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

        iDonor = 0;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iFace = 0; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {

            faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
            iNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1]- (unsigned int)faceindex;
            for (iTarget=0; iTarget<iNodes; iTarget++) {
              if (Global_Point == Buffer_Receive_GlobalPoint[faceindex+iTarget]) {
                coeff =Buffer_Receive_Coeff[faceindex+iTarget];
                pGlobalPoint = Buffer_Receive_FaceNodes[faceindex+iTarget];
                target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint(iDonor,pGlobalPoint);
                target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(iDonor,coeff);
                target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, iProcessor);
                //cout <<rank << " Global Point " << Global_Point<<" iDonor " << iDonor <<" coeff " << coeff <<" gp " << pGlobalPoint << endl;
                iDonor++;
              }
            }
          }
        }
      }
    }
    delete[] Buffer_Send_nFace_Donor;
    delete[] Buffer_Send_nFaceNodes_Donor;

    delete[] Buffer_Receive_nFace_Donor;
    delete[] Buffer_Receive_nFaceNodes_Donor;

    delete[] Buffer_Send_FaceIndex;
    delete[] Buffer_Send_FaceNodes;
    delete[] Buffer_Send_GlobalPoint;
    delete[] Buffer_Send_Coeff;

    delete[] Buffer_Receive_FaceIndex;
    delete[] Buffer_Receive_FaceNodes;
    delete[] Buffer_Receive_GlobalPoint;
    delete[] Buffer_Receive_Coeff;

  }
}

CSlidingMesh::CSlidingMesh(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone)  :  CInterpolator(geometry_container, config, iZone, jZone){

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

  /*--- For fluid-structure interaction data interpolated with have nDim dimensions ---*/
 // InitializeData(Zones,nDim);
}

CSlidingMesh::~CSlidingMesh(){}

void CSlidingMesh::Set_TransferCoeff(CConfig **config){
    
  /* --- This routine sets the transfer coefficient for sliding mesh approach --- */
  
  /*
   * The algorithm is based on Rinaldi et al. "Flux-conserving treatment of non-conformal interfaces 
   * for finite-volume discritization of conservaation laws" 2015, Comp. Fluids, 120, pp 126-139
   */

  /*  0 - Variable declaration - */

  /* --- General variables --- */

  bool check;
  
  unsigned short iDim, nDim;
  
  unsigned long ii, jj, *uptr;
  unsigned long vPoint, dPoint;
  unsigned long iEdgeVisited, nEdgeVisited, iNodeVisited;
  unsigned long nAlreadyVisited, nToVisit, StartVisited;
  
  unsigned long *alreadyVisitedDonor, *ToVisit, *tmpVect;
  unsigned long *storeProc, *tmp_storeProc;

  su2double dTMP;
  su2double *Coeff_Vect, *tmp_Coeff_Vect;               

  /* --- Geometrical variables --- */

  su2double *Coord_i, *Coord_j, dist, mindist, *Normal;
  su2double Area, Area_old, tmp_Area;
  su2double LineIntersectionLength, *Direction, length;


  /* --- Markers Variables --- */

  unsigned short iMarkerInt, nMarkerInt; 

  unsigned long iVertex, nVertexTarget;

  int markDonor, markTarget;

  /* --- Target variables --- */

  unsigned long target_iPoint, jVertexTarget;
  unsigned long nEdges_target, nNode_target;

  unsigned long *Target_nLinkedNodes, *Target_LinkedNodes, *Target_StartLinkedNodes, *target_segment;
  unsigned long *Target_GlobalPoint, *Target_Proc;  
  
  su2double *TargetPoint_Coord, *target_iMidEdge_point, *target_jMidEdge_point, **target_element;

  /* --- Donor variables --- */

  unsigned long donor_StartIndex, donor_forward_point, donor_backward_point, donor_iPoint, donor_OldiPoint; 
  unsigned long nEdges_donor, nNode_donor, nGlobalVertex_Donor; 

  unsigned long nDonorPoints, iDonor;
  unsigned long *Donor_Vect, *tmp_Donor_Vect;
  unsigned long *Donor_nLinkedNodes, *Donor_LinkedNodes, *Donor_StartLinkedNodes;
  unsigned long *Donor_GlobalPoint, *Donor_Proc;
  
  su2double *donor_iMidEdge_point, *donor_jMidEdge_point;
  su2double **donor_element, *DonorPoint_Coord;
    
  /*  1 - Variable pre-processing - */

  nDim = donor_geometry->GetnDim();

  /*--- Setting up auxiliary vectors ---*/

  Donor_Vect = NULL;
  Coeff_Vect = NULL;
  storeProc  = NULL;

  tmp_Donor_Vect = NULL;
  tmp_Coeff_Vect = NULL;
  tmp_storeProc  = NULL;
  
  Normal    = new su2double[nDim];
  Direction = new su2double[nDim];
    
    
  /* 2 - Find boundary tag between touching grids */

  /*--- Number of markers on the FSI interface ---*/
  nMarkerInt    = (int)( config[ donorZone ]->GetMarker_n_ZoneInterface() ) / 2;

  /*--- For the number of markers on the interface... ---*/
  for ( iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++ ){

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if( !CheckInterfaceBoundary(markDonor, markTarget) )
      continue;

    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;

    /*
    3 -Reconstruct the boundaries from parallel partitioning
    */

    /*--- Target boundary ---*/
    ReconstructBoundary(targetZone, markTarget);
    
    nGlobalVertex_Target = nGlobalVertex;

    TargetPoint_Coord       = Buffer_Receive_Coord;
    Target_GlobalPoint      = Buffer_Receive_GlobalPoint;
    Target_nLinkedNodes     = Buffer_Receive_nLinkedNodes;
    Target_StartLinkedNodes = Buffer_Receive_StartLinkedNodes;
    Target_LinkedNodes      = Buffer_Receive_LinkedNodes;
    Target_Proc             = Buffer_Receive_Proc;
    
    /*--- Donor boundary ---*/
    ReconstructBoundary(donorZone, markDonor);
    
    nGlobalVertex_Donor = nGlobalVertex;

    DonorPoint_Coord       = Buffer_Receive_Coord;
    Donor_GlobalPoint      = Buffer_Receive_GlobalPoint;
    Donor_nLinkedNodes     = Buffer_Receive_nLinkedNodes;
    Donor_StartLinkedNodes = Buffer_Receive_StartLinkedNodes;
    Donor_LinkedNodes      = Buffer_Receive_LinkedNodes;
    Donor_Proc             = Buffer_Receive_Proc;

    /*--- Starts building the supermesh layer (2D or 3D) ---*/
    /* - For each target node, it first finds the closest donor point
     * - Then it creates the supermesh in the close proximity of the target point:
     * - Starting from the closest donor node, it expands the supermesh by including 
     * donor elements neighboring the initial one, until the overall target area is fully covered.
     */

    if(nDim == 2){
        
      target_iMidEdge_point = new su2double[nDim];
      target_jMidEdge_point = new su2double[nDim];

      donor_iMidEdge_point = new su2double[nDim];
      donor_jMidEdge_point = new su2double[nDim];
        
      /*--- Starts with supermesh reconstruction ---*/

      target_segment = new unsigned long[2];
      
      for (iVertex = 0; iVertex < nVertexTarget; iVertex++) {

        nDonorPoints = 0;

        /*--- Stores coordinates of the target node ---*/

        target_iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();

        if (target_geometry->node[target_iPoint]->GetDomain()){

          Coord_i = target_geometry->node[target_iPoint]->GetCoord();

          /*--- Brute force to find the closest donor_node ---*/

          mindist = 1E6;
          donor_StartIndex = 0;
 
          for (donor_iPoint = 0; donor_iPoint < nGlobalVertex_Donor; donor_iPoint++) {
        
            Coord_j = &DonorPoint_Coord[ donor_iPoint * nDim ];

            dist = PointsDistance(Coord_i, Coord_j);

            if (dist < mindist) {
              mindist = dist;  
              donor_StartIndex = donor_iPoint;
            }

            if (dist == 0.0){
              donor_StartIndex = donor_iPoint;
              break;
            }    
          }

          donor_iPoint    = donor_StartIndex;
          donor_OldiPoint = donor_iPoint;
          
          /*--- Contruct information regarding the target cell ---*/
          
          dPoint = target_geometry->node[target_iPoint]->GetGlobalIndex();
          for (jVertexTarget = 0; jVertexTarget < nGlobalVertex_Target; jVertexTarget++)
            if( dPoint == Target_GlobalPoint[jVertexTarget] )
              break;
            
          if ( Target_nLinkedNodes[jVertexTarget] == 1 ){
            target_segment[0] = Target_LinkedNodes[ Target_StartLinkedNodes[jVertexTarget] ];
            target_segment[1] = jVertexTarget;
          }
          else{
            target_segment[0] = Target_LinkedNodes[ Target_StartLinkedNodes[jVertexTarget] ];
            target_segment[1] = Target_LinkedNodes[ Target_StartLinkedNodes[jVertexTarget] + 1];
          }
      
          dTMP = 0;
          for(iDim = 0; iDim < nDim; iDim++){
            target_iMidEdge_point[iDim] = ( TargetPoint_Coord[ nDim * target_segment[0] + iDim ] + target_geometry->node[ target_iPoint ]->GetCoord(iDim) ) / 2;
            target_jMidEdge_point[iDim] = ( TargetPoint_Coord[ nDim * target_segment[1] + iDim ] + target_geometry->node[ target_iPoint ]->GetCoord(iDim) ) / 2;

            Direction[iDim] = target_jMidEdge_point[iDim] - target_iMidEdge_point[iDim];
            dTMP += Direction[iDim] * Direction[iDim];
          }

          dTMP = sqrt(dTMP);
          for(iDim = 0; iDim < nDim; iDim++)
            Direction[iDim] /= dTMP;

          length = PointsDistance(target_iMidEdge_point, target_jMidEdge_point);

          check = false;

          /*--- Proceeds along the forward direction (depending on which connected boundary node is found first) ---*/

          while( !check ){
  
            /*--- Proceeds until the value of the intersection area is null ---*/

            if ( Donor_nLinkedNodes[donor_iPoint] == 1 ){
              donor_forward_point  = Donor_LinkedNodes[ Donor_StartLinkedNodes[donor_iPoint] ];
              donor_backward_point = donor_iPoint;
            }
            else{
              uptr = &Donor_LinkedNodes[ Donor_StartLinkedNodes[donor_iPoint] ];
              
              if( donor_OldiPoint != uptr[0] ){
                donor_forward_point  = uptr[0];
                donor_backward_point = uptr[1];
              }
              else{
                donor_forward_point  = uptr[1];
                donor_backward_point = uptr[0];
              }
            }
            
            if(donor_iPoint >= nGlobalVertex_Donor){
              check = true;
              continue;
            }
            
            for(iDim = 0; iDim < nDim; iDim++){
              donor_iMidEdge_point[iDim] = ( DonorPoint_Coord[ donor_forward_point  * nDim + iDim] + DonorPoint_Coord[ donor_iPoint * nDim + iDim] ) / 2;
              donor_jMidEdge_point[iDim] = ( DonorPoint_Coord[ donor_backward_point * nDim + iDim] + DonorPoint_Coord[ donor_iPoint * nDim + iDim] ) / 2;
            }

            LineIntersectionLength = ComputeLineIntersectionLength(target_iMidEdge_point, target_jMidEdge_point, donor_iMidEdge_point, donor_jMidEdge_point, Direction);

            if ( LineIntersectionLength == 0.0 ){
              check = true;
              continue;
            }
  
            /*--- In case the element intersects the target cell, update the auxiliary communication data structure ---*/

            tmp_Coeff_Vect = new     su2double[ nDonorPoints + 1 ];
            tmp_Donor_Vect = new unsigned long[ nDonorPoints + 1 ];
            tmp_storeProc  = new unsigned long[ nDonorPoints + 1 ];
 
            for( iDonor = 0; iDonor < nDonorPoints; iDonor++){
              tmp_Donor_Vect[iDonor] = Donor_Vect[iDonor];
              tmp_Coeff_Vect[iDonor] = Coeff_Vect[iDonor];
              tmp_storeProc[iDonor]  = storeProc[iDonor];
            }

            tmp_Donor_Vect[ nDonorPoints ] = donor_iPoint;
            tmp_Coeff_Vect[ nDonorPoints ] = LineIntersectionLength / length;
            tmp_storeProc[  nDonorPoints ] = Donor_Proc[donor_iPoint];
            
            if (Donor_Vect != NULL) delete [] Donor_Vect;  
            if (Coeff_Vect != NULL) delete [] Coeff_Vect;            
            if (storeProc  != NULL) delete [] storeProc;

            Donor_Vect = tmp_Donor_Vect;
            Coeff_Vect = tmp_Coeff_Vect;
            storeProc  = tmp_storeProc;

            donor_OldiPoint = donor_iPoint;
            donor_iPoint    = donor_forward_point;

            nDonorPoints++;
          }
             
          if ( Donor_nLinkedNodes[donor_StartIndex] == 2 ){
            check = false;
           
            uptr = &Donor_LinkedNodes[ Donor_StartLinkedNodes[donor_StartIndex] ];

            donor_iPoint = uptr[1];
            donor_OldiPoint = donor_StartIndex;
          }
          else
            check = true;

          /*--- Proceeds along the backward direction (depending on which connected boundary node is found first) ---*/

          while( !check ){

            /*--- Proceeds until the value of the intersection length is null ---*/
            if ( Donor_nLinkedNodes[donor_iPoint] == 1 ){
              donor_forward_point  = donor_OldiPoint;
              donor_backward_point = donor_iPoint;
            }
            else{
              uptr = &Donor_LinkedNodes[ Donor_StartLinkedNodes[donor_iPoint] ];
              
              if( donor_OldiPoint != uptr[0] ){
                donor_forward_point  = uptr[0];
                donor_backward_point = uptr[1];
              }
              else{
                donor_forward_point  = uptr[1];
                donor_backward_point = uptr[0];
              }
            }

            if(donor_iPoint >= nGlobalVertex_Donor){
              check = true;
              continue;
            }
            
            for(iDim = 0; iDim < nDim; iDim++){
              donor_iMidEdge_point[iDim] = ( DonorPoint_Coord[ donor_forward_point  * nDim + iDim] + DonorPoint_Coord[ donor_iPoint * nDim + iDim] ) / 2;
              donor_jMidEdge_point[iDim] = ( DonorPoint_Coord[ donor_backward_point * nDim + iDim] + DonorPoint_Coord[ donor_iPoint * nDim + iDim] ) / 2;
            }       

            LineIntersectionLength = ComputeLineIntersectionLength(target_iMidEdge_point, target_jMidEdge_point, donor_iMidEdge_point, donor_jMidEdge_point, Direction);

            if ( LineIntersectionLength == 0.0 ){
              check = true;
              continue;
            }

            /*--- In case the element intersects the target cell, update the auxiliary communication data structure ---*/

            tmp_Coeff_Vect = new     su2double[ nDonorPoints + 1 ];
            tmp_Donor_Vect = new unsigned long[ nDonorPoints + 1 ];
            tmp_storeProc  = new unsigned long[ nDonorPoints + 1 ];
 
            for( iDonor = 0; iDonor < nDonorPoints; iDonor++){
              tmp_Donor_Vect[iDonor] = Donor_Vect[iDonor];
              tmp_Coeff_Vect[iDonor] = Coeff_Vect[iDonor];
              tmp_storeProc[iDonor]  = storeProc[iDonor];
            }
                  
            tmp_Coeff_Vect[ nDonorPoints ] = LineIntersectionLength / length;                  
            tmp_Donor_Vect[ nDonorPoints ] = donor_iPoint;
            tmp_storeProc[  nDonorPoints ] = Donor_Proc[donor_iPoint];

            if (Donor_Vect != NULL) delete [] Donor_Vect;
            if (Coeff_Vect != NULL) delete [] Coeff_Vect;
            if (storeProc  != NULL) delete [] storeProc;

            Donor_Vect = tmp_Donor_Vect;
            Coeff_Vect = tmp_Coeff_Vect;
            storeProc  = tmp_storeProc;

            donor_OldiPoint = donor_iPoint;
            donor_iPoint    = donor_forward_point;
  
            nDonorPoints++;
          }
                
          /*--- Set the communication data structure and copy data from the auxiliary vectors ---*/

          target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nDonorPoints);

          target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();
   
          for ( iDonor = 0; iDonor < nDonorPoints; iDonor++ ){              
            target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(          iDonor, Coeff_Vect[iDonor]);
            target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint(    iDonor, Donor_GlobalPoint[ Donor_Vect[iDonor] ]);
            target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, storeProc[iDonor]);
          }
        }
      }    
      
      delete [] target_segment;
      
      delete [] target_iMidEdge_point;
      delete [] target_jMidEdge_point;

      delete [] donor_iMidEdge_point;
      delete [] donor_jMidEdge_point;
    }
    else{ 
      /* --- 3D geometry, creates a superficial super-mesh --- */
      
      for (iVertex = 0; iVertex < nVertexTarget; iVertex++) {
  
        nDonorPoints = 0;

        /*--- Stores coordinates of the target node ---*/

        target_iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();
        
        if (target_geometry->node[target_iPoint]->GetDomain()){
    
          Coord_i = target_geometry->node[target_iPoint]->GetCoord();

          target_geometry->vertex[markTarget][iVertex]->GetNormal(Normal);
 
          /*--- The value of Area computed here includes also portion of boundary belonging to different marker ---*/
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) 
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
 
          for (iDim = 0; iDim < nDim; iDim++)
            Normal[iDim] /= Area;

          for (iDim = 0; iDim < nDim; iDim++)
            Coord_i[iDim] = target_geometry->node[target_iPoint]->GetCoord(iDim);
          
          dPoint = target_geometry->node[target_iPoint]->GetGlobalIndex();
          for (target_iPoint = 0; target_iPoint < nGlobalVertex_Target; target_iPoint++){
            if( dPoint == Target_GlobalPoint[target_iPoint] )
              break;
          }        
        
          /*--- Build local surface dual mesh for target element ---*/
        
          nEdges_target = Target_nLinkedNodes[target_iPoint];

          nNode_target = 2*(nEdges_target + 1);
          
          target_element = new su2double*[nNode_target];
          for (ii = 0; ii < nNode_target; ii++)
            target_element[ii] = new su2double[nDim];
            
          nNode_target = Build_3D_surface_element(Target_LinkedNodes, Target_StartLinkedNodes, Target_nLinkedNodes, TargetPoint_Coord, target_iPoint, target_element);

          /*--- Brute force to find the closest donor_node ---*/

          mindist = 1E6;
          donor_StartIndex = 0;
 
          for (donor_iPoint = 0; donor_iPoint < nGlobalVertex_Donor; donor_iPoint++) {
        
            Coord_j = &DonorPoint_Coord[ donor_iPoint * nDim ];

            dist = PointsDistance(Coord_i, Coord_j);

            if (dist < mindist) {
              mindist = dist;  
              donor_StartIndex = donor_iPoint;
            }

            if (dist == 0.0){
              donor_StartIndex = donor_iPoint;
              break;
            }    
          }
                
          donor_iPoint = donor_StartIndex;

          nEdges_donor = Donor_nLinkedNodes[donor_iPoint];

          donor_element = new su2double*[ 2*nEdges_donor + 2 ];
          for (ii = 0; ii < 2*nEdges_donor + 2; ii++)
            donor_element[ii] = new su2double[nDim];                

          nNode_donor = Build_3D_surface_element(Donor_LinkedNodes, Donor_StartLinkedNodes, Donor_nLinkedNodes, DonorPoint_Coord, donor_iPoint, donor_element);

          Area = 0;
          for (ii = 1; ii < nNode_target-1; ii++){
            for (jj = 1; jj < nNode_donor-1; jj++){
              Area += Compute_Triangle_Intersection(target_element[0], target_element[ii], target_element[ii+1], donor_element[0], donor_element[jj], donor_element[jj+1], Normal);
              //cout << Compute_Triangle_Intersection(target_element[0], target_element[ii], target_element[ii+1], donor_element[0], donor_element[jj], donor_element[jj+1], Normal) << endl;
            }
          }

          for (ii = 0; ii < 2*nEdges_donor + 2; ii++)
            delete [] donor_element[ii];
          delete [] donor_element;

          nDonorPoints = 1;

          /*--- In case the element intersect the target cell update the auxiliary communication data structure ---*/

          Coeff_Vect = new     su2double[ nDonorPoints ];
          Donor_Vect = new unsigned long[ nDonorPoints ];
          storeProc  = new unsigned long[ nDonorPoints ];

          Coeff_Vect[0] = Area;
          Donor_Vect[0] = donor_iPoint;
          storeProc[0]  = Donor_Proc[donor_iPoint];

          alreadyVisitedDonor = new unsigned long[1];

          alreadyVisitedDonor[0] = donor_iPoint;
          nAlreadyVisited = 1;
          StartVisited = 0;

          Area_old = -1;
                    
          while( Area > Area_old ){ 

            /* 
             * - Starting from the closest donor_point, it expands the supermesh by a countour search pattern.
             * - The closest donor element becomes the core, at each iteration a new layer of elements around the core is taken into account
             */

            Area_old = Area;

            ToVisit = NULL;
            nToVisit = 0;

            for( iNodeVisited = StartVisited; iNodeVisited < nAlreadyVisited; iNodeVisited++ ){

              vPoint = alreadyVisitedDonor[ iNodeVisited ];
            
              nEdgeVisited = Donor_nLinkedNodes[vPoint];
 
              for (iEdgeVisited = 0; iEdgeVisited < nEdgeVisited; iEdgeVisited++){

                donor_iPoint = Donor_LinkedNodes[ Donor_StartLinkedNodes[vPoint] + iEdgeVisited];

                /*--- Check if the node to visit is already listed in the data structure to avoid double visits ---*/

                check = 0;

                for( jj = 0; jj < nAlreadyVisited; jj++ ){
                  if( donor_iPoint == alreadyVisitedDonor[jj] ){
                    check = 1; 
                    break;
                  }
                }

                if( check == 0 && ToVisit != NULL){
                  for( jj = 0; jj < nToVisit; jj++ )
                    if( donor_iPoint == ToVisit[jj] ){
                      check = 1; 
                      break;
                    }       
                }

                if( check == 0 ){ 
                  /*--- If the node was not already visited, visit it and list it into data structure ---*/
  
                  tmpVect = new unsigned long[ nToVisit + 1 ];

                  for( jj = 0; jj < nToVisit; jj++ )
                    tmpVect[jj] = ToVisit[jj];
                  tmpVect[nToVisit] = donor_iPoint;

                  if( ToVisit != NULL )
                    delete [] ToVisit;
                    
                  ToVisit = tmpVect;
                  tmpVect = NULL;

                  nToVisit++; 

                  /*--- Find the value of the intersection area between the current donor element and the target element --- */

                  nEdges_donor = Donor_nLinkedNodes[donor_iPoint];

                  donor_element = new su2double*[ 2*nEdges_donor + 2 ];   
                  for (ii = 0; ii < 2*nEdges_donor + 2; ii++)
                    donor_element[ii] = new su2double[nDim];             

                  nNode_donor = Build_3D_surface_element(Donor_LinkedNodes, Donor_StartLinkedNodes, Donor_nLinkedNodes, DonorPoint_Coord, donor_iPoint, donor_element);

                  tmp_Area = 0;
                  for (ii = 1; ii < nNode_target-1; ii++)
                    for (jj = 1; jj < nNode_donor-1; jj++)
                      tmp_Area += Compute_Triangle_Intersection(target_element[0], target_element[ii], target_element[ii+1], donor_element[0], donor_element[jj], donor_element[jj+1], Normal);

                  for (ii = 0; ii < 2*nEdges_donor + 2; ii++)
                    delete [] donor_element[ii];
                  delete [] donor_element;
 
                  /*--- In case the element intersect the target cell update the auxiliary communication data structure ---*/

                  tmp_Coeff_Vect = new     su2double[ nDonorPoints + 1 ];
                  tmp_Donor_Vect = new unsigned long[ nDonorPoints + 1 ];
                  tmp_storeProc  = new unsigned long[ nDonorPoints + 1 ];
 
                  for( iDonor = 0; iDonor < nDonorPoints; iDonor++){
                    tmp_Donor_Vect[iDonor] = Donor_Vect[iDonor];
                    tmp_Coeff_Vect[iDonor] = Coeff_Vect[iDonor];
                    tmp_storeProc[iDonor]  = storeProc[iDonor];
                  }
                  
                  tmp_Coeff_Vect[ nDonorPoints ] = tmp_Area;                  
                  tmp_Donor_Vect[ nDonorPoints ] = donor_iPoint;
                  tmp_storeProc[  nDonorPoints ] = Donor_Proc[donor_iPoint];

                  if (Donor_Vect != NULL) {delete [] Donor_Vect; }
                  if (Coeff_Vect != NULL) {delete [] Coeff_Vect; }
                  if (storeProc  != NULL) {delete [] storeProc;  }

                  Donor_Vect = tmp_Donor_Vect;
                  Coeff_Vect = tmp_Coeff_Vect;
                  storeProc  = tmp_storeProc;

                  tmp_Coeff_Vect = NULL;                  
                  tmp_Donor_Vect = NULL;
                  tmp_storeProc  = NULL;

                  nDonorPoints++;
   
                  Area += tmp_Area;
                }
              }   
            }

            /*--- Update auxiliary data structure ---*/
 
            StartVisited = nAlreadyVisited;

            tmpVect = new unsigned long[ nAlreadyVisited + nToVisit ];

            for( jj = 0; jj < nAlreadyVisited; jj++ )
              tmpVect[jj] = alreadyVisitedDonor[jj];
              
            for( jj = 0; jj < nToVisit; jj++ )
              tmpVect[ nAlreadyVisited + jj ] = ToVisit[jj];

            if( alreadyVisitedDonor != NULL )
              delete [] alreadyVisitedDonor;

            alreadyVisitedDonor = tmpVect;            

            nAlreadyVisited += nToVisit;    

            delete [] ToVisit;  
          }

          delete [] alreadyVisitedDonor;
  
          /*--- Set the communication data structure and copy data from the auxiliary vectors ---*/

          target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nDonorPoints);
          target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

          for ( iDonor = 0; iDonor < nDonorPoints; iDonor++ ){              
            target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(iDonor, Coeff_Vect[iDonor]/Area);
            target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint( iDonor, Donor_GlobalPoint[ Donor_Vect[iDonor] ] );
            target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, storeProc[iDonor]);
            //cout <<rank << " Global Point " << Global_Point<<" iDonor " << iDonor <<" coeff " << coeff <<" gp " << pGlobalPoint << endl;               
          }

          for (ii = 0; ii < 2*nEdges_target + 2; ii++)
            delete [] target_element[ii];
          delete [] target_element;
          
          if (Donor_Vect != NULL) {delete [] Donor_Vect; Donor_Vect = NULL;}
          if (Coeff_Vect != NULL) {delete [] Coeff_Vect; Coeff_Vect = NULL;}
          if (storeProc  != NULL) {delete [] storeProc;  storeProc  = NULL;}
        }
      }
    }


    delete [] TargetPoint_Coord;
    delete [] Target_GlobalPoint;
    delete [] Target_Proc;
    delete [] Target_nLinkedNodes;
    delete [] Target_LinkedNodes;
    delete [] Target_StartLinkedNodes;
        
    delete [] DonorPoint_Coord;
    delete [] Donor_GlobalPoint;
    delete [] Donor_Proc;
    delete [] Donor_nLinkedNodes;      
    delete [] Donor_StartLinkedNodes;  
    delete [] Donor_LinkedNodes;       
    
  }

  delete [] Normal;
  delete [] Direction;
  
  if (Donor_Vect != NULL) delete [] Donor_Vect;
  if (Coeff_Vect != NULL) delete [] Coeff_Vect;
  if (storeProc  != NULL) delete [] storeProc;  
}

int CSlidingMesh::Build_3D_surface_element(unsigned long *map, unsigned long *startIndex, unsigned long* nNeighbor, su2double *coord, unsigned long centralNode, su2double** element){
    
  /*--- Given a node "centralNode", this routines reconstruct the vertex centered surface element around the node and store it into "element" ---*/
  /*--- Returns the number of points included in the element ---*/
  
  unsigned long iNode, jNode, kNode, iElementNode, iPoint, jPoint, nOuterNodes;

  unsigned short nDim = 3, iDim, nTmp;

  int NextNode, **OuterNodesNeighbour, CurrentNode, StartIndex, count;
  unsigned long *OuterNodes, *ptr;

  /* --- Store central node as element first point --- */

  for (iDim = 0; iDim < nDim; iDim++)
    element[0][iDim] = coord[centralNode * nDim + iDim];

  nOuterNodes = nNeighbor[centralNode];

  OuterNodes = &map[ startIndex[centralNode] ];

  /* --- Allocate auxiliary structure, vectors are longer than needed but this avoid further re-allocations due to length variation --- */

  OuterNodesNeighbour = new int*[nOuterNodes];
  for ( iNode = 0; iNode < nOuterNodes; iNode++ )
    OuterNodesNeighbour[ iNode ] = new int[2];

  /* --- Finds which and how many nodes belong to the specified marker, initialize some variables --- */

  for ( iNode = 0; iNode < nOuterNodes; iNode++ ){
    OuterNodesNeighbour[ iNode ][0] = -1;
    OuterNodesNeighbour[ iNode ][1] = -1;
  }
    
  /* --- For each outer node, the program finds the two neighbouring outer nodes --- */

  StartIndex = 0;
  for( iNode = 0; iNode < nOuterNodes; iNode++ ){

  count = 0;
  iPoint = OuterNodes[ iNode ];
  ptr = &map[ startIndex[iPoint] ];
  nTmp = nNeighbor[iPoint];

  for ( jNode = 0; jNode < nTmp; jNode++ ){
    jPoint = ptr[jNode];
    for( kNode = 0; kNode < nOuterNodes; kNode++ ){
      if ( jPoint == OuterNodes[ kNode ] && jPoint != centralNode){
        OuterNodesNeighbour[iNode][count] = (int)kNode;
        count++;
        break;
      }
    }
  }

  // If the central node belongs to two different markers, ie at corners, makes this outer node the starting point for reconstructing the element
  if( count == 1 ) 
    StartIndex = (int)iNode;
  }

  /* --- Build element, starts from one outer node and loops along the external edges until the element is reconstructed --- */

  CurrentNode = StartIndex;
  NextNode    = OuterNodesNeighbour[ CurrentNode ][0];
  iElementNode = 1;

  while( NextNode != -1 ){

    for (iDim = 0; iDim < nDim; iDim++)
      element[ iElementNode ][iDim] = ( element[0][iDim] + coord[ OuterNodes[ CurrentNode ] * nDim + iDim ])/2;

    iElementNode++;

    for (iDim = 0; iDim < nDim; iDim++) 
      element[ iElementNode ][iDim] = ( element[0][iDim] + coord[ OuterNodes[ CurrentNode ] * nDim + iDim] + coord[ OuterNodes[ NextNode ] * nDim + iDim] )/3;

    iElementNode++;

    if( OuterNodesNeighbour[ NextNode ][0] == CurrentNode){
      CurrentNode = NextNode; 
      NextNode = OuterNodesNeighbour[ NextNode ][1];  
    }
    else{
      CurrentNode = NextNode; 
      NextNode = OuterNodesNeighbour[ NextNode ][0];  
    }

    if (CurrentNode == StartIndex)
      break;
    }

    if( CurrentNode == StartIndex ){ // This is a closed element, so add again element 1 to the end of the structure, useful later

    for (iDim = 0; iDim < nDim; iDim++)
      element[ iElementNode ][iDim] = element[1][iDim];
    iElementNode++;
  }
  else{
    for (iDim = 0; iDim < nDim; iDim++)
    element[ iElementNode ][iDim] = ( element[0][iDim] + coord[ OuterNodes[ CurrentNode ] * nDim + iDim] )/2;
    iElementNode++;
  }

  for ( iNode = 0; iNode < nOuterNodes; iNode++ )
    delete [] OuterNodesNeighbour[ iNode ];
  delete [] OuterNodesNeighbour;

  return (int)iElementNode;
  
}

su2double CSlidingMesh::ComputeLineIntersectionLength(su2double* A1, su2double* A2, su2double* B1, su2double* B2, su2double* Direction){
    
  /*--- Given 2 segments, each defined by 2 points, it projects them along a given direction and it computes the length of the segment resulting from their intersection ---*/
  /*--- The algorithm works for both 2D and 3D problems ---*/

  unsigned short iDim;
  unsigned short nDim = donor_geometry->GetnDim();

  su2double dotA2, dotB1, dotB2;

  dotA2 = 0;
  for(iDim = 0; iDim < nDim; iDim++)
    dotA2 += ( A2[iDim] - A1[iDim] ) * Direction[iDim];

  if( dotA2 >= 0 ){
    dotB1 = 0;
    dotB2 = 0;
    for(iDim = 0; iDim < nDim; iDim++){
      dotB1 += ( B1[iDim] - A1[iDim] ) * Direction[iDim];
      dotB2 += ( B2[iDim] - A1[iDim] ) * Direction[iDim];
    }
  }
  else{
    dotA2 *= -1;

    dotB1 = 0;
    dotB2 = 0;
    for(iDim = 0; iDim < nDim; iDim++){
      dotB1 -= ( B1[iDim] - A1[iDim] ) * Direction[iDim];
      dotB2 -= ( B2[iDim] - A1[iDim] ) * Direction[iDim];
    }
  }

  if( dotB1 >= 0 && dotB1 <= dotA2 ){
    if ( dotB2 < 0 )
      return fabs( dotB1 );
    if ( dotB2 > dotA2 )
      return fabs( dotA2 - dotB1 );

    return fabs( dotB1 - dotB2 );
  }

  if( dotB2 >= 0 && dotB2 <= dotA2 ){
    if ( dotB1 < 0 )
      return fabs(dotB2);
    if ( dotB1 > dotA2 )
      return fabs( dotA2 - dotB2 );
  }

  if( ( dotB1 <= 0 && dotA2 <= dotB2 ) || ( dotB2 <= 0 && dotA2 <= dotB1 ) )
    return fabs( dotA2 );

  return 0.0;
}

su2double CSlidingMesh::Compute_Triangle_Intersection(su2double* A1, su2double* A2, su2double* A3, su2double* B1, su2double* B2, su2double* B3, su2double* Direction){
    
  /* --- This routine is ONLY for 3D grids --- */
  /* --- Projects triangle points onto a plane, specified by its normal "Direction", and calls the ComputeIntersectionArea routine --- */

  unsigned short iDim;
  unsigned short nDim = 3;

  su2double I[3], J[3], K[3];
  su2double a1[3], a2[3], a3[3];
  su2double b1[3], b2[3], b3[3];
  su2double m1, m2;

  /* --- Reference frame is determined by: x = A1A2 y = x ^ ( -Direction ) --- */

  for(iDim = 0; iDim < 3; iDim++){
    a1[iDim] = 0;
    a2[iDim] = 0;
    a3[iDim] = 0;

    b1[iDim] = 0;
    b2[iDim] = 0;
    b3[iDim] = 0;
  }

  m1 = 0;
  for(iDim = 0; iDim < nDim; iDim++){
    K[iDim] = Direction[iDim];

    m1 += K[iDim] * K[iDim];
  }

  for(iDim = 0; iDim < nDim; iDim++)
    K[iDim] /= sqrt(m1);

  m2 = 0;
  for(iDim = 0; iDim < nDim; iDim++)
    m2 += (A2[iDim] - A1[iDim]) * K[iDim];

  m1 = 0;
  for(iDim = 0; iDim < nDim; iDim++){
    I[iDim] = (A2[iDim] - A1[iDim]) - m2 * K[iDim];
    m1 += I[iDim] * I[iDim];
  }

  for(iDim = 0; iDim < nDim; iDim++)
    I[iDim] /= sqrt(m1);

  // Cross product to find Y
  J[0] =   K[1]*I[2] - K[2]*I[1];
  J[1] = -(K[0]*I[2] - K[2]*I[0]);
  J[2] =   K[0]*I[1] - K[1]*I[0];

  /* --- Project all points on the plane specified by Direction and change their reference frame taking A1 as origin --- */

  for(iDim = 0; iDim < nDim; iDim++){
    a2[0] += (A2[iDim] - A1[iDim]) * I[iDim];
    a2[1] += (A2[iDim] - A1[iDim]) * J[iDim];
    a2[2] += (A2[iDim] - A1[iDim]) * K[iDim];

    a3[0] += (A3[iDim] - A1[iDim]) * I[iDim];
    a3[1] += (A3[iDim] - A1[iDim]) * J[iDim];
    a3[2] += (A3[iDim] - A1[iDim]) * K[iDim];

    b1[0] += (B1[iDim] - A1[iDim]) * I[iDim];
    b1[1] += (B1[iDim] - A1[iDim]) * J[iDim];
    b1[2] += (B1[iDim] - A1[iDim]) * K[iDim];

    b2[0] += (B2[iDim] - A1[iDim]) * I[iDim];
    b2[1] += (B2[iDim] - A1[iDim]) * J[iDim];
    b2[2] += (B2[iDim] - A1[iDim]) * K[iDim];

    b3[0] += (B3[iDim] - A1[iDim]) * I[iDim];
    b3[1] += (B3[iDim] - A1[iDim]) * J[iDim];
    b3[2] += (B3[iDim] - A1[iDim]) * K[iDim];       
  }

  /* --- Compute intersection area --- */

  return ComputeIntersectionArea( a1, a2, a3, b1, b2, b3 );
}

su2double CSlidingMesh::ComputeIntersectionArea( su2double* P1, su2double* P2, su2double* P3, su2double* Q1, su2double* Q2, su2double* Q3 ){
    
  /* --- This routines computes the area of the polygonal element generated by the superimposition of 2 planar triangle --- */
  /* --- The 2 triangle must lie on the same plane --- */

  unsigned short iDim, nPoints, i, j, k;
  unsigned short nDim, min_theta_index;

  su2double points[16][2], IntersectionPoint[2], theta[6];
  su2double TriangleP[4][2], TriangleQ[4][2];
  su2double Area, det, dot1, dot2, dtmp, min_theta;

  nDim    = 2;
  nPoints = 0;

  for(iDim = 0; iDim < nDim; iDim++){
    TriangleP[0][iDim] = 0;
    TriangleP[1][iDim] = P2[iDim] - P1[iDim];
    TriangleP[2][iDim] = P3[iDim] - P1[iDim];
    TriangleP[3][iDim] = 0;

    TriangleQ[0][iDim] = Q1[iDim] - P1[iDim];
    TriangleQ[1][iDim] = Q2[iDim] - P1[iDim];
    TriangleQ[2][iDim] = Q3[iDim] - P1[iDim];
    TriangleQ[3][iDim] = Q1[iDim] - P1[iDim];
  }


  for( j = 0; j < 3; j++){
    if( CheckPointInsideTriangle(TriangleP[j], TriangleQ[0], TriangleQ[1], TriangleQ[2]) ){

      // Then P1 is also inside triangle Q, so store it
      for(iDim = 0; iDim < nDim; iDim++)
        points[nPoints][iDim] = TriangleP[j][iDim];

      nPoints++;      
    }
  }

  for( j = 0; j < 3; j++){    
    if( CheckPointInsideTriangle(TriangleQ[j], TriangleP[0], TriangleP[1], TriangleP[2]) ){

      // Then Q1 is also inside triangle P, so store it
      for(iDim = 0; iDim < nDim; iDim++)
        points[nPoints][iDim] = TriangleQ[j][iDim];

      nPoints++;      
    }
  }


  // Compute all edge intersections

  for( j = 0; j < 3; j++){
    for( i = 0; i < 3; i++){

      det = (TriangleP[j][0] - TriangleP[j+1][0]) * ( TriangleQ[i][1] - TriangleQ[i+1][1] ) - (TriangleP[j][1] - TriangleP[j+1][1]) * (TriangleQ[i][0] - TriangleQ[i+1][0]);

      if ( det != 0.0 ){
        ComputeLineIntersectionPoint( TriangleP[j], TriangleP[j+1], TriangleQ[i], TriangleQ[i+1], IntersectionPoint );

        dot1 = 0;
        dot2 = 0;
        for(iDim = 0; iDim < nDim; iDim++){
          dot1 += ( TriangleP[j][iDim] - IntersectionPoint[iDim] ) * ( TriangleP[j+1][iDim] - IntersectionPoint[iDim] );
          dot2 += ( TriangleQ[i][iDim] - IntersectionPoint[iDim] ) * ( TriangleQ[i+1][iDim] - IntersectionPoint[iDim] );
        }

       if( dot1 <= 0 && dot2 <= 0 ){ // It found one intersection

         // Store temporarily the intersection point

         for(iDim = 0; iDim < nDim; iDim++)
           points[nPoints][iDim] = IntersectionPoint[iDim];

         nPoints++;
       }   
     }
   }
 }

  // Remove double points, if any

  for( i = 0; i < nPoints; i++){
    for( j = i+1; j < nPoints; j++){
      if(points[j][0] == points[i][0] && points[j][1] == points[i][1]){
        for( k = j; k < nPoints-1; k++){
          points[k][0] = points[k+1][0];
          points[k][1] = points[k+1][1];
        }
        nPoints--;
        j--;
      }
    }
  }       

  // Re-order nodes   

  for( i = 1; i < nPoints; i++){ // Change again reference frame
    for(iDim = 0; iDim < nDim; iDim++)
      points[i][iDim] -= points[0][iDim]; 

    // Compute polar azimuth for each node but the first
    theta[i] = atan2(points[i][1], points[i][0]);
  }

  for(iDim = 0; iDim < nDim; iDim++)
    points[0][iDim] = 0;

  for( i = 1; i < nPoints; i++){

    min_theta = theta[i];
    min_theta_index = 0;

    for( j = i + 1; j < nPoints; j++){ 

      if( theta[j] < min_theta ){
        min_theta = theta[j];
        min_theta_index = j;
      }
    }

    if( min_theta_index != 0 ){
      dtmp = theta[i];
      theta[i] = theta[min_theta_index];
      theta[min_theta_index] = dtmp;

      dtmp = points[i][0];
      points[i][0] = points[min_theta_index][0];
      points[min_theta_index][0] = dtmp;

      dtmp = points[i][1];
      points[i][1] = points[min_theta_index][1];
      points[min_theta_index][1] = dtmp;
    }
  }
    
  // compute area using cross product rule, points position are referred to the 2-dimensional, local, reference frame centered in points[0]

  Area = 0;

  if (nPoints > 2){
    for( i = 1; i < nPoints-1; i++ ){

      // Ax*By
      Area += ( points[i][0] - points[0][0] ) * ( points[i+1][1] - points[0][1] );

      // Ay*Bx
      Area -= ( points[i][1] - points[0][1] ) * ( points[i+1][0] - points[0][0] );        
    }
  }

  return fabs(Area)/2;
}

void CSlidingMesh::ComputeLineIntersectionPoint( su2double* A1, su2double* A2, su2double* B1, su2double* B2, su2double* IntersectionPoint ){
    
  /* --- Uses determinant rule to compute the intersection point between 2 straight segments --- */
  /* This works only for lines on a 2D plane, A1, A2 and B1, B2 are respectively the head and the tail points of each segment, 
   * since they're on a 2D plane they are defined by a 2-elements array containing their coordinates */

  su2double det;

  det = (A1[0] - A2[0]) * (B1[1] - B2[1]) - (A1[1] - A2[1]) * (B1[0] - B2[0]);
 
  if ( det != 0.0 ){ // else there is no intersection point
    IntersectionPoint[0] = ( ( A1[0]*A2[1] - A1[1]*A2[0] ) * ( B1[0] - B2[0] ) - ( B1[0]*B2[1] - B1[1]*B2[0] ) * ( A1[0] - A2[0] ) ) / det;
    IntersectionPoint[1] = ( ( A1[0]*A2[1] - A1[1]*A2[0] ) * ( B1[1] - B2[1] ) - ( B1[0]*B2[1] - B1[1]*B2[0] ) * ( A1[1] - A2[1] ) ) / det;
  }

  return;
}

bool CSlidingMesh::CheckPointInsideTriangle(su2double* Point, su2double* T1, su2double* T2, su2double* T3){

  /* --- Check whether a point "Point" lies inside or outside a triangle defined by 3 points "T1", "T2", "T3" --- */
  /* For each edge it checks on which side the point lies:
   * - Computes the unit vector pointing at the internal side of the edge
   * - Comutes the vector that connects the point to a point along the edge
   * - If the dot product is positive it means that the point is on the internal side of the edge
   * - If the check is positive for all the 3 edges, then the point lies within the triangle
   */ 

  unsigned short iDim, nDim, check;

  su2double vect1[2], vect2[2], r[2];
  su2double dot;

  check = 0;
  nDim  = 2;

  /* --- Check first edge --- */

  dot = 0;
  for(iDim = 0; iDim < nDim; iDim++){
    vect1[iDim] = T3[iDim] - T1[iDim]; // vec 1 is aligned to the edge
    vect2[iDim] = T2[iDim] - T1[iDim]; // vect 2 is the vector connecting one edge point to the third triangle vertex

    r[iDim] = Point[iDim] - T1[iDim];  // Connects point to vertex T1

    dot += vect2[iDim] * vect2[iDim];  
  }
  dot = sqrt(dot);

  for(iDim = 0; iDim < nDim; iDim++)
    vect2[iDim] /= dot;

  dot = 0;
  for(iDim = 0; iDim < nDim; iDim++)
    dot += vect1[iDim] * vect2[iDim];

  for(iDim = 0; iDim < nDim; iDim++)
    vect1[iDim] = T3[iDim] - (T1[iDim] + dot * vect2[iDim]); // Computes the inward unit vector

  dot = 0;
  for(iDim = 0; iDim < nDim; iDim++)  // Checs that the point lies on the internal plane
    dot += vect1[iDim] * r[iDim];

  if (dot >= 0)
    check++;

  /* --- Check second edge --- */

  dot = 0;
  for(iDim = 0; iDim < nDim; iDim++){
    vect1[iDim] = T1[iDim] - T2[iDim];
    vect2[iDim] = T3[iDim] - T2[iDim];

    r[iDim] = Point[iDim] - T2[iDim];

    dot += vect2[iDim] * vect2[iDim];
  }
  dot = sqrt(dot);

  for(iDim = 0; iDim < nDim; iDim++)
    vect2[iDim] /= dot;

  dot = 0;
  for(iDim = 0; iDim < nDim; iDim++)
    dot += vect1[iDim] * vect2[iDim];

  for(iDim = 0; iDim < nDim; iDim++)
    vect1[iDim] = T1[iDim] - (T2[iDim] + dot * vect2[iDim]);

  dot = 0;
  for(iDim = 0; iDim < nDim; iDim++)
    dot += vect1[iDim] * r[iDim];

  if (dot >= 0)
    check++;

  /* --- Check third edge --- */

  dot = 0;
  for(iDim = 0; iDim < nDim; iDim++){
    vect1[iDim] = T2[iDim] - T3[iDim];
    vect2[iDim] = T1[iDim] - T3[iDim];

    r[iDim] = Point[iDim] - T3[iDim];

    dot += vect2[iDim] * vect2[iDim];
  }
  dot = sqrt(dot);

  for(iDim = 0; iDim < nDim; iDim++)
    vect2[iDim] /= dot;

  dot = 0;
  for(iDim = 0; iDim < nDim; iDim++)
    dot += vect1[iDim] * vect2[iDim];

  for(iDim = 0; iDim < nDim; iDim++)
    vect1[iDim] = T2[iDim] - (T3[iDim] + dot * vect2[iDim]);

  dot = 0;
  for(iDim = 0; iDim < nDim; iDim++)
    dot += vect1[iDim] * r[iDim];

  if (dot >= 0)
    check++;

  return (check == 3);     
}
