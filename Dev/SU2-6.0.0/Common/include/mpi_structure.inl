/*!
 * \file mpi_structure.hpp
 * \brief In-Line subroutines of the <i>mpi_structure.hpp</i> file.
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

#pragma once

#ifdef HAVE_MPI

inline void CBaseMPIWrapper::Error(std::string ErrorMsg, std::string FunctionName){
  if (Rank == 0){
    std::cout << std::endl << std::endl;
    std::cout << "Error in \"" << FunctionName << "\": " << std::endl;
    std::cout <<  "-------------------------------------------------------------------------" << std::endl;
    std::cout << ErrorMsg << std::endl;
    std::cout <<  "------------------------------ Error Exit -------------------------------" << std::endl;
    std::cout << std::endl << std::endl;    
  }
  Abort(currentComm, 0);
}


inline int CBaseMPIWrapper::GetRank(){
  return Rank;
}

inline int CBaseMPIWrapper::GetSize(){
  return Size;
}

inline void CBaseMPIWrapper::SetComm(Comm newComm){
  currentComm = newComm;
  MPI_Comm_rank(currentComm, &Rank);  
  MPI_Comm_size(currentComm, &Size);
}

inline CBaseMPIWrapper::Comm CBaseMPIWrapper::GetComm(){
  return currentComm;
}

inline void CBaseMPIWrapper::Init(int *argc, char ***argv) {
  MPI_Init(argc,argv);
  MPI_Comm_rank(currentComm, &Rank);    
  MPI_Comm_size(currentComm, &Size);  
}

inline void CBaseMPIWrapper::Buffer_attach(void *buffer, int size){
  MPI_Buffer_attach(buffer, size);
}

inline void CBaseMPIWrapper::Buffer_detach(void *buffer, int *size){
  MPI_Buffer_detach(buffer, size);
}

inline void CBaseMPIWrapper::Comm_rank(Comm comm, int *rank){
  MPI_Comm_rank(comm, rank);
}

inline void CBaseMPIWrapper::Comm_size(Comm comm, int *size){
  MPI_Comm_size(comm, size);
}

inline void CBaseMPIWrapper::Finalize(){
  MPI_Finalize();
}

inline void CBaseMPIWrapper::Barrier(Comm comm) {
  MPI_Barrier(comm);
}

inline void CBaseMPIWrapper::Abort(Comm comm, int error) {
  MPI_Abort(comm, error);
}

inline void CBaseMPIWrapper::Get_count(Status *status, Datatype datatype, int *count) {
  MPI_Get_count(status, datatype, count);
}

inline void CBaseMPIWrapper::Isend(void *buf, int count, Datatype datatype,
                               int dest, int tag, Comm comm, Request *request) {
  MPI_Isend(buf,count,datatype,dest,tag,comm,request);
}

inline void CBaseMPIWrapper::Irecv(void *buf, int count, Datatype datatype,
                               int dest, int tag, Comm comm, Request *request) {
  MPI_Irecv(buf,count,datatype,dest,tag,comm, request);
}

inline void CBaseMPIWrapper::Wait(Request *request, Status *status) {
  MPI_Wait(request,status);
}

inline void CBaseMPIWrapper::Testall(int count, Request *array_of_requests, int *flag, Status *array_of_statuses) {
  MPI_Testall(count,array_of_requests,flag, array_of_statuses);
}

inline void CBaseMPIWrapper::Waitall(int nrequests, Request *request, Status *status) {
  MPI_Waitall(nrequests, request, status);
}

inline void CBaseMPIWrapper::Probe(int source, int tag, Comm comm, Status *status){
  MPI_Probe(source, tag, comm, status);
}

inline void CBaseMPIWrapper::Send(void *buf, int count, Datatype datatype,
                              int dest, int tag, Comm comm) {
  MPI_Send(buf,count,datatype,dest,tag,comm);
}

inline void CBaseMPIWrapper::Recv(void *buf, int count, Datatype datatype,
                              int dest,int tag, Comm comm, Status *status) {
  MPI_Recv(buf,count,datatype,dest,tag,comm,status);
}

inline void CBaseMPIWrapper::Bcast(void *buf, int count, Datatype datatype,
                               int root, Comm comm) {
  MPI_Bcast(buf,count,datatype,root,comm);
}

inline void CBaseMPIWrapper::Bsend(void *buf, int count, Datatype datatype,
                               int dest, int tag, Comm comm) {
  MPI_Bsend(buf,count,datatype,dest,tag,comm);
}

inline void CBaseMPIWrapper::Reduce(void *sendbuf, void *recvbuf, int count,
                                Datatype datatype, Op op, int root, Comm comm) {
  MPI_Reduce(sendbuf, recvbuf,count,datatype,op,root,comm);
}

inline void CBaseMPIWrapper::Allreduce(void *sendbuf, void *recvbuf, int count,
                                   Datatype datatype, Op op, Comm comm) {
  MPI_Allreduce(sendbuf,recvbuf,count,datatype,op,comm);
}

inline void CBaseMPIWrapper::Gather(void *sendbuf, int sendcnt,Datatype sendtype,
                                void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm) {
  MPI_Gather(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,root,comm);
}

inline void CBaseMPIWrapper::Scatter(void *sendbuf, int sendcnt,Datatype sendtype,
                                 void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm) {
  MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
}

inline void CBaseMPIWrapper::Allgather(void *sendbuf, int sendcnt, Datatype sendtype,
                                   void *recvbuf, int recvcnt, Datatype recvtype, Comm comm) {
  MPI_Allgather(sendbuf,sendcnt,sendtype, recvbuf, recvcnt, recvtype, comm);
}

inline void CBaseMPIWrapper::Allgatherv(void *sendbuf, int sendcount, Datatype sendtype,
                                        void *recvbuf, int *recvcounts, int *displs, Datatype recvtype, Comm comm){
  MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
}

inline void CBaseMPIWrapper::Alltoall(void *sendbuf, int sendcount, Datatype sendtype, void *recvbuf, int recvcount, Datatype recvtype, Comm comm){
  MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}

inline void CBaseMPIWrapper::Sendrecv(void *sendbuf, int sendcnt, Datatype sendtype,
                                  int dest, int sendtag, void *recvbuf, int recvcnt,
                                  Datatype recvtype,int source, int recvtag,
                                  Comm comm, Status *status) {
  MPI_Sendrecv(sendbuf,sendcnt,sendtype,dest,sendtag,recvbuf,recvcnt,recvtype,source,recvtag,comm,status);
}

inline void CBaseMPIWrapper::Waitany(int nrequests, Request *request,
                                 int *index, Status *status) {
  MPI_Waitany(nrequests, request, index, status);
}
  

#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE

inline void CMediMPIWrapper::Init(int *argc, char ***argv) {
  AMPI_Init(argc,argv);
  MediTool::init();
  AMPI_Comm_rank(convertComm(currentComm), &Rank);    
  AMPI_Comm_size(convertComm(currentComm), &Size);  
}

inline void CMediMPIWrapper::SetComm(Comm newComm){
  currentComm = newComm;
  AMPI_Comm_rank(convertComm(currentComm), &Rank);  
  AMPI_Comm_size(convertComm(currentComm), &Size);
}

inline AMPI_Comm CMediMPIWrapper::convertComm(MPI_Comm comm) {
  return comm;
}

inline AMPI_Datatype CMediMPIWrapper::convertDatatype(MPI_Datatype datatype){
  if (datatype == MPI_DOUBLE){
    return AMPI_ADOUBLE;
  }
  else if (datatype == MPI_SHORT){
    return AMPI_SHORT;
  }
  else if (datatype == MPI_UNSIGNED_SHORT){
    return AMPI_UNSIGNED_SHORT;
  }
  else if (datatype == MPI_LONG){
    return AMPI_LONG;
  }
  else if (datatype == MPI_UNSIGNED_LONG){
    return AMPI_UNSIGNED_LONG;
  }
  else if (datatype == MPI_SHORT){
    return AMPI_SHORT;
  }
  else if (datatype == MPI_CHAR){
    return AMPI_CHAR;
  }
  else if (datatype == MPI_INT){
    return AMPI_INT;
  } else {
    Error("Conversion not implemented", CURRENT_FUNCTION);
    return AMPI_DOUBLE;
  }
}

inline AMPI_Op CMediMPIWrapper::convertOp(MPI_Op op) {
  if(MPI_SUM == op) {
    return medi::AMPI_SUM;
  } else if(MPI_PROD == op) {
    return medi::AMPI_PROD;
  } else if(MPI_MIN == op) {
    return medi::AMPI_MIN;
  } else if(MPI_MAX == op) {
    return medi::AMPI_MAX;
  } else {
    Error("Conversion not implemented", CURRENT_FUNCTION);
    return medi::AMPI_SUM;
  }
}
inline void CMediMPIWrapper::Buffer_attach(void *buffer, int size){
  AMPI_Buffer_attach(buffer, size);
}

inline void CMediMPIWrapper::Buffer_detach(void *buffer, int *size){
  AMPI_Buffer_detach(buffer, size);
}

inline void CMediMPIWrapper::Comm_rank(Comm comm, int *rank){
  AMPI_Comm_rank(convertComm(comm), rank);
}

inline void CMediMPIWrapper::Comm_size(Comm comm, int *size){
  AMPI_Comm_size(convertComm(comm), size);
}

inline void CMediMPIWrapper::Finalize(){
  AMPI_Finalize();
}

inline void CMediMPIWrapper::Barrier(Comm comm) {
  AMPI_Barrier(convertComm(comm));
}

inline void CMediMPIWrapper::Abort(Comm comm, int error) {
  AMPI_Abort(convertComm(comm), error);
}

inline void CMediMPIWrapper::Get_count(Status *status, Datatype datatype, int *count) {
  AMPI_Get_count(status, convertDatatype(datatype), count);
}

inline void CMediMPIWrapper::Isend(void *buf, int count, Datatype datatype,
                               int dest, int tag, Comm comm, Request *request) {
  AMPI_Isend(buf,count,convertDatatype(datatype),dest,tag,convertComm(comm),request);
}

inline void CMediMPIWrapper::Irecv(void *buf, int count, Datatype datatype,
                               int dest, int tag, Comm comm, Request *request) {
  AMPI_Irecv(buf,count,convertDatatype(datatype),dest,tag,convertComm(comm), request);
}

inline void CMediMPIWrapper::Wait(SU2_MPI::Request *request, Status *status) {
  AMPI_Wait(request,status);
}

inline void CMediMPIWrapper::Testall(int count, Request *array_of_requests, int *flag, Status *array_of_statuses) {
  AMPI_Testall(count,array_of_requests,flag, array_of_statuses);
}

inline void CMediMPIWrapper::Waitall(int nrequests, Request *request, Status *status) {
  AMPI_Waitall(nrequests, request, status);
}

inline void CMediMPIWrapper::Probe(int source, int tag, Comm comm, Status *status){
  AMPI_Probe(source, tag, convertComm(comm), status);
}

inline void CMediMPIWrapper::Send(void *buf, int count, Datatype datatype,
                              int dest, int tag, Comm comm) {
  AMPI_Send(buf,count,convertDatatype(datatype),dest,tag,convertComm(comm));
}

inline void CMediMPIWrapper::Recv(void *buf, int count, Datatype datatype,
                              int dest,int tag, Comm comm, Status *status) {
  AMPI_Recv(buf,count,convertDatatype(datatype),dest,tag,convertComm(comm),status);
}

inline void CMediMPIWrapper::Bcast(void *buf, int count, Datatype datatype,
                               int root, Comm comm) {
  AMPI_Bcast(buf,count,convertDatatype(datatype),root,convertComm(comm));
}

inline void CMediMPIWrapper::Bsend(void *buf, int count, Datatype datatype,
                               int dest, int tag, Comm comm) {
  AMPI_Bsend(buf,count,convertDatatype(datatype),dest,tag,convertComm(comm));
}

inline void CMediMPIWrapper::Reduce(void *sendbuf, void *recvbuf, int count,
                                Datatype datatype, Op op, int root, Comm comm) {
  AMPI_Reduce(sendbuf, recvbuf,count,convertDatatype(datatype),convertOp(op),root,convertComm(comm));
}

inline void CMediMPIWrapper::Allreduce(void *sendbuf, void *recvbuf, int count,
                                   Datatype datatype, Op op, Comm comm) {
  AMPI_Allreduce(sendbuf,recvbuf,count,convertDatatype(datatype),convertOp(op),convertComm(comm));
}

inline void CMediMPIWrapper::Gather(void *sendbuf, int sendcnt,Datatype sendtype,
                                void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm) {
  AMPI_Gather(sendbuf,sendcnt,convertDatatype(sendtype),recvbuf,recvcnt,convertDatatype(recvtype),root,convertComm(comm));
}

inline void CMediMPIWrapper::Scatter(void *sendbuf, int sendcnt,Datatype sendtype,
                                 void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm) {
  AMPI_Scatter(sendbuf, sendcnt, convertDatatype(sendtype), recvbuf, recvcnt, convertDatatype(recvtype), root, convertComm(comm));
}

inline void CMediMPIWrapper::Allgather(void *sendbuf, int sendcnt, Datatype sendtype,
                                   void *recvbuf, int recvcnt, Datatype recvtype, Comm comm) {
  AMPI_Allgather(sendbuf,sendcnt,convertDatatype(sendtype), recvbuf, recvcnt, convertDatatype(recvtype), convertComm(comm));
}

inline void CMediMPIWrapper::Allgatherv(void *sendbuf, int sendcount, Datatype sendtype,
                                        void *recvbuf, int *recvcounts, int *displs, Datatype recvtype, Comm comm){
  AMPI_Allgatherv(sendbuf, sendcount, convertDatatype(sendtype), recvbuf, recvcounts, displs, convertDatatype(recvtype), convertComm(comm));
}

inline void CMediMPIWrapper::Alltoall(void *sendbuf, int sendcount, Datatype sendtype, void *recvbuf, int recvcount, Datatype recvtype, Comm comm){
  AMPI_Alltoall(sendbuf, sendcount, convertDatatype(sendtype), recvbuf, recvcount, convertDatatype(recvtype), convertComm(comm));
}

inline void CMediMPIWrapper::Sendrecv(void *sendbuf, int sendcnt, Datatype sendtype,
                                  int dest, int sendtag, void *recvbuf, int recvcnt,
                                  Datatype recvtype,int source, int recvtag,
                                  Comm comm, Status *status) {
  AMPI_Sendrecv(sendbuf,sendcnt,convertDatatype(sendtype),dest,sendtag,recvbuf,recvcnt,convertDatatype(recvtype),source,recvtag,convertComm(comm),status);
}

inline void CMediMPIWrapper::Waitany(int nrequests, Request *request,
                                 int *index, Status *status) {
  AMPI_Waitany(nrequests, request, index, status);
}
#endif
#else

inline void CBaseMPIWrapper::Error(std::string ErrorMsg, std::string FunctionName){
  if (Rank == 0){
    std::cout << std::endl << std::endl;
    std::cout << "Error in \"" << FunctionName << "\": " << std::endl;
    std::cout <<  "-------------------------------------------------------------------------" << std::endl;
    std::cout << ErrorMsg << std::endl;
    std::cout <<  "------------------------------ Error Exit -------------------------------" << std::endl;
    std::cout << std::endl << std::endl;    
  }
  Abort(currentComm, 0);
}

inline int CBaseMPIWrapper::GetRank(){
  return Rank;
}

inline int CBaseMPIWrapper::GetSize(){
  return Size;
}

inline void CBaseMPIWrapper::SetComm(Comm newComm){
  currentComm = newComm;
}

inline CBaseMPIWrapper::Comm CBaseMPIWrapper::GetComm(){
  return currentComm;
}

inline void CBaseMPIWrapper::Init(int *argc, char ***argv) {}

inline void CBaseMPIWrapper::Buffer_attach(void *buffer, int size) {}

inline void CBaseMPIWrapper::Buffer_detach(void *buffer, int *size) {}

inline void CBaseMPIWrapper::Barrier(Comm comm) {}

inline void CBaseMPIWrapper::Abort(Comm comm, int error) {exit(EXIT_FAILURE);}

inline void CBaseMPIWrapper::Comm_rank(Comm comm, int *rank) {*rank = 0;}

inline void CBaseMPIWrapper::Comm_size(Comm comm, int *size) {*size = 1;}

inline void CBaseMPIWrapper::Finalize(){}

inline void CBaseMPIWrapper::Isend(void *buf, int count, Datatype datatype, int dest,
                               int tag, Comm comm, Request* request) {}

inline void CBaseMPIWrapper::Irecv(void *buf, int count, Datatype datatype, int source,
                               int tag, Comm comm, Request* request) {}

inline void CBaseMPIWrapper::Wait(Request *request, Status *status) {}

inline void CBaseMPIWrapper::Waitall(int nrequests, Request *request, Status *status) {}

inline void CBaseMPIWrapper::Waitany(int nrequests, Request *request,
                                 int *index, Status *status) {}

inline void CBaseMPIWrapper::Send(void *buf, int count, Datatype datatype, int dest,
                              int tag, Comm comm) {}

inline void CBaseMPIWrapper::Recv(void *buf, int count, Datatype datatype, int dest,
                              int tag, Comm comm, Status *status) {}

inline void CBaseMPIWrapper::Bcast(void *buf, int count, Datatype datatype, int root,
                               Comm comm) {}

inline void CBaseMPIWrapper::Bsend(void *buf, int count, Datatype datatype, int dest,
                               int tag, Comm comm){}

inline void  CBaseMPIWrapper::Reduce(void *sendbuf, void *recvbuf, int count,
                                 Datatype datatype, Op op, int root, Comm comm){
   CopyData(sendbuf, recvbuf, count, datatype);
}

inline void  CBaseMPIWrapper::Allreduce(void *sendbuf, void *recvbuf, int count,
                                    Datatype datatype, Op op, Comm comm){
   CopyData(sendbuf, recvbuf, count, datatype);
}

inline void CBaseMPIWrapper::Gather(void *sendbuf, int sendcnt, Datatype sendtype,
                   void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm){
  CopyData(sendbuf, recvbuf, sendcnt, sendtype);

}

inline void CBaseMPIWrapper::Scatter(void *sendbuf, int sendcnt, Datatype sendtype,
                    void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm){
  CopyData(sendbuf, recvbuf, sendcnt, sendtype);

}

inline void CBaseMPIWrapper::Allgatherv(void *sendbuf, int sendcnt, Datatype sendtype,
                                   void *recvbuf, int recvcnt, int *displs, Datatype recvtype, Comm comm){
  CopyData(sendbuf, recvbuf, sendcnt, sendtype);
}

inline void CBaseMPIWrapper::Allgather(void *sendbuf, int sendcnt, Datatype sendtype,
                                   void *recvbuf, int recvcnt, Datatype recvtype, Comm comm){
  CopyData(sendbuf, recvbuf, sendcnt, sendtype);

}

inline void CBaseMPIWrapper::Sendrecv(void *sendbuf, int sendcnt, Datatype sendtype,
                                  int dest, int sendtag, void *recvbuf, int recvcnt,
                                  Datatype recvtype,int source, int recvtag,
                                  Comm comm, Status *status){
  CopyData(sendbuf, recvbuf, sendcnt, sendtype);

}

inline void CBaseMPIWrapper::Alltoall(void *sendbuf, int sendcount, Datatype sendtype,
                                  void *recvbuf, int recvcount, Datatype recvtype,
                                  Comm comm){
  CopyData(sendbuf, recvbuf, recvcount, sendtype);
}

inline void CBaseMPIWrapper::Probe(int source, int tag, Comm comm, Status *status){}

inline void CBaseMPIWrapper::CopyData(void *sendbuf, void *recvbuf, int size, Datatype datatype){
  switch (datatype) {
    case MPI_DOUBLE:
      for (int i = 0; i < size; i++){
        static_cast<su2double*>(recvbuf)[i] = static_cast<su2double*>(sendbuf)[i];
      }
      break;
    case MPI_UNSIGNED_LONG:
      for (int i = 0; i < size; i++){
        static_cast<unsigned long*>(recvbuf)[i] = static_cast<unsigned long*>(sendbuf)[i];
      }
      break;
    case MPI_LONG:
      for (int i = 0; i < size; i++){
        static_cast<long*>(recvbuf)[i] = static_cast<long*>(sendbuf)[i];
      }
      break;
    case MPI_UNSIGNED_SHORT:
      for (int i = 0; i < size; i++){
        static_cast<unsigned short*>(recvbuf)[i] = static_cast<unsigned short*>(sendbuf)[i];
      }
      break;
    case MPI_CHAR:
      for (int i = 0; i < size; i++){
        static_cast<char*>(recvbuf)[i] = static_cast<char*>(sendbuf)[i];
      }
      break;
    case MPI_SHORT:
      for (int i = 0; i < size; i++){
        static_cast<short*>(recvbuf)[i] = static_cast<short*>(sendbuf)[i];
      }
      break;
    case MPI_INT:
      for (int i = 0; i < size; i++){
        static_cast<int*>(recvbuf)[i] = static_cast<int*>(sendbuf)[i];
      }
      break;
    default:
      break;
  }
}
#endif
