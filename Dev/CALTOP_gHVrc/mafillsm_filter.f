!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine mafillsm_filter(co,kon,ipkon,lakon,ne,ttime,time,
     &  mortar,ne0,nea,neb,elCentroid)
!
!     filling the stiffness matrix in spare matrix format (sm)
!
      implicit none
!
c      integer 
!
      character*8 lakon(*)
!
      integer kon(*),ne0,
     &  ipkon(*),ne,mortar,
     &  nea,neb,i
!
      real*8 co(3,*),ttime,
     &  time,xcg,ycg,zcg,elCentroid(3,*)
!
      intent(in) co,kon,ipkon,lakon,ne,
     &  ttime,time,
     &  mortar,ne0,nea,neb
!
      intent(inout) elCentroid

c      penall=3.0
c      write(*,*) loc(kflag)
c      write(*,*) loc(s)
c      write(*,*) loc(sm)
c      write(*,*) loc(ff)
c      write(*,*) loc(index1)
!

!
c       open (unit=200,file="Mafillsmfilter.dat")
 
!
!
!     mechanical analysis: loop over all elements
!
      do i=nea,neb
           call e_c3d_filter(co,kon,lakon(i),
     &          i,
     &          ttime,time,ne0,ipkon,mortar,
     &          xcg,ycg,zcg)

            elCentroid(1,i)=xcg
            elCentroid(2,i)=ycg
            elCentroid(3,i)=zcg
c      write(200,*) i,xcg,ycg,zcg
      enddo

 

c     close(200)
!
      return
      end
