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
      subroutine e_c3d_filter2(co,kon,lakonl,
     &  nelem,
     &  ttime,time,ne0,ipkon,mortar,xcg,ycg,zcg)

!
!     computation of the element centroid for the element with
!     the topology in konl
!
!
      implicit none
!
      character*8 lakonl
!
      integer konl(26),nelem,
     &  mortar,kon(*),i,j,k,l,kl,ipkon(*),indexe,
     &  nope,nopes,nopev,ne0
!
      real*8 co(3,*),xl(3,26),
     &  xl2(3,9),coords(3),xl1(3,9),dtime,ttime,time,
     &  tvar(2),xcg,ycg,zcg
!
      intent(in) co,kon,lakonl,
     &  nelem,ttime,time,ne0,ipkon,mortar
!
      intent(inout) xcg,ycg,zcg

!
      tvar(1)=time
      tvar(2)=ttime+time
      xcg=0.d0
      ycg=0.d0
      zcg=0.d0
!
      open (unit=300,file="e_c3dfilter.dat")
      indexe=ipkon(nelem)
c     Bernhardi start
      if(lakonl(1:5).eq.'C3D8I') then
         nope=11
         nopev=8
         nopes=4
      elseif(lakonl(4:5).eq.'20') then
c     Bernhardi end
         nope=20
         nopev=8
         nopes=8
      elseif(lakonl(4:4).eq.'8') then
         nope=8
         nopev=8
         nopes=4
      elseif(lakonl(4:5).eq.'10') then
         nope=10
         nopev=4
         nopes=6
      elseif(lakonl(4:4).eq.'4') then
         nope=4
         nopev=4
         nopes=3
      elseif(lakonl(4:5).eq.'15') then
         nope=15
         nopev=6
      elseif(lakonl(4:4).eq.'6') then
         nope=6
         nopev=6
      elseif(lakonl(1:2).eq.'ES') then
         if(lakonl(7:7).eq.'C') then
            if(mortar.eq.0) then
               nope=ichar(lakonl(8:8))-47
               konl(nope+1)=kon(indexe+nope+1)
            elseif(mortar.eq.1) then
               nope=kon(indexe)
            endif
         else
            nope=ichar(lakonl(8:8))-47
         endif
      elseif(lakonl(1:4).eq.'MASS') then
         nope=1
      endif
!

      do i=1,nope
        konl(i)=kon(indexe+i)
        do j=1,3
          xl(j,i)=co(j,konl(i))
          write(300,*) i,j,xl(j,i)
        enddo
!      calculate cg based on xl coordinates
          xcg=xcg+xl(1,i)
          ycg=ycg+xl(2,i)
          zcg=zcg+xl(3,i)      
      enddo

       xcg=xcg/nope
       ycg=ycg/nope
       zcg=zcg/nope


!
!  
      close(300)
      return
      end

