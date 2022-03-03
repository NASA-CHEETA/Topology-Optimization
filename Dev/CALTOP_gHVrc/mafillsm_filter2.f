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
      subroutine mafillsm_filter2(ne,ttime,time,
     &  ne0,nea,neb,elCentroid,
     &  rmin,filternnz,
     &  FilterMatrixs,rowFilters,colFilters,filternnzElems,elarr,
     &  fnnzassumed)
!
!     calculating the distance between element centroids
!
      implicit none
!
!
!
      integer ne0,filternnz,
     &  ne,nea,neb,i,j,elarr(*),ii, fnnzassumed,
     &  filternnzElems(*),rowFilters(fnnzassumed,*),
     &  colFilters(fnnzassumed,*), dummy1
!
      real*8 FilterMatrixs(fnnzassumed,*),ttime,
     &  time,rmin,elCentroid(3,*),sum,
     &  xi,yi,zi,xj,yj,zj,d,rmind
!
      intent(in) ne,elCentroid,rmin,
     &  ttime,time,fnnzassumed,
     &  ne0,nea,neb,elarr
!
      intent(inout) FilterMatrixs,
     &  filternnz,
     &  rowFilters,colFilters,filternnzElems

c      penall=3.0
c      write(*,*) loc(kflag)
c      write(*,*) loc(s)
c      write(*,*) loc(sm)
c      write(*,*) loc(ff)
c      write(*,*) loc(index1)
!
	dummy1=fnnzassumed/3
!
c      open (unit=200,file="Mafillsmfilter2.dat")
 
!
!       
        
!      loop over all elements
!
cc      do i=nea,neb
!      centroid of element i
cc       xi=elCentroid(1,i)
cc       yi=elCentroid(2,i)
cc       zi=elCentroid(3,i)

cc      sum=0
cc       filternnz=0
cc       filternnzElem(i)=0
cc        do j=1,ne   
!        centroid of element j
cc         xj=elCentroid(1,j)
cc         yj=elCentroid(2,j)
cc         zj=elCentroid(3,j)
cc          d=((xj-xi)**2+(yj-yi)**2+(zj-zi)**2)**0.5
cc          rmind=rmin-d

cc          if(rmind.GT.0) then
cc            filternnz=filternnz+1
cc            filternnzElem(i)=filternnzElem(i)+1

cc            rowFilter(filternnz,i)=i
cc            colFilter(filternnz,i)=j
cc            FilterMatrix(filternnz,i)=rmind**2
cc            sum=sum+rmind**2
cc          endif

cc        enddo
c      Divide by sum
cc      do j=1,filternnz
cc      FilterMatrix(j,i)=FilterMatrix(j,i)/sum
                                
cc      enddo

c      write(200,*) i,xcg,ycg,zcg
cc      enddo

c       write(200,*) FilterMatrix,"...",colFilter,rowFilter

c      close(200)
!


       do ii=nea,neb
c      get element number, indexed from 0 in C 
       i=elarr(ii)+1
!      centroid of element i
       xi=elCentroid(1,i)
       yi=elCentroid(2,i)
       zi=elCentroid(3,i)

       sum=0
       filternnz=0
       filternnzElems(i)=0
         do j=i,ne   
!         centroid of element j
          xj=elCentroid(1,j)
          yj=elCentroid(2,j)
          zj=elCentroid(3,j)
           d=((xj-xi)**2+(yj-yi)**2+(zj-zi)**2)
           rmind=rmin*rmin-d

             if(rmind.GE.0) then
               filternnz=filternnz+1
               filternnzElems(i)=filternnzElems(i)+1

                rowFilters(filternnz,i)=i
                colFilters(filternnz,i)=j
                FilterMatrixs(filternnz,i)=rmin-d**0.5

              endif

	     if(filternnzElems(i).EQ.(dummy1)) then
		EXIT
	     endif

           enddo

      enddo

      return
      end
