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
      subroutine mafillsmvectorfilter(ne,ttime,time,
     &  ne0,nea,neb,FilterMatrix,Vector,VectorFiltered,
     &  filternnzElem,rowFilter,colFilter,fnnzassumed,q)
!
!     filtering Vector array to VectorFiltered
!
      implicit none
!
c      integer 
!
!
      integer filternnzElem(*),fnnzassumed, ne0,row,col,
     &  ne,nea,neb,i,j,rowFilter(fnnzassumed,*),
     &  colFilter(fnnzassumed,*)
!
      real*8 FilterMatrix(fnnzassumed,*),ttime,sum,
     &  time, Vector(*),VectorFiltered(*),q
 
!
      intent(in) ne,FilterMatrix, fnnzassumed,
     &  ttime,time,filternnzElem,Vector,
     &  ne0,nea,neb,rowFilter,colFilter,q
!
      intent(inout) VectorFiltered

c      penall=3.0
c      write(*,*) loc(kflag)
c      write(*,*) loc(s)
c      write(*,*) loc(sm)
c      write(*,*) loc(ff)
c      write(*,*) loc(index1)
!

!
c       open (unit=500,file="MMafillsm_Vectorfilter.dat")
 
!
!       
        
!      loop over all elements
!

      
      do i=nea,neb
        sum=0.d0    
        do j=1,filternnzElem(i)  
        VectorFiltered(i)=VectorFiltered(i)+
     &   (FilterMatrix(j,i)**q)*Vector(ColFilter(j,i))

            sum=sum+FilterMatrix(j,i)**q

 
        enddo
c       write(500,*) i, Vector(i), VectorFiltered(i)
        VectorFiltered(i)=VectorFiltered(i)/sum    
      enddo

c       write(500,*) FilterMatrix,"...",colFilter,rowFilter

c      close(500)
!
      return
      end
