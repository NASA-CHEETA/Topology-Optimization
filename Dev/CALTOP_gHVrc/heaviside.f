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
      subroutine heaviside(ne,nea, neb, designFiltered,designFilteredH,
     &  dHSVector,
     &  beta_H, eta_H)
!
!     filtering with Heaviside projection
!
      implicit none
!
      integer ne,
     &  nea,neb,i
!
      real*8 designFiltered(*),dHSVector(*),beta_H,eta_H,x,
     &  designFilteredH(*), fmin
 
!
      intent(in) ne,nea,neb,beta_H,eta_H, designFiltered

      intent(inout) designFilteredH, dHSVector

!      loop over all elements
!
      fmin = 0.d0
      do i=nea,neb
        x = designFiltered(i)
        designFilteredH(i)=(tanh(beta_H*eta_H)
     &            +tanh(beta_H*(x-eta_H)))
     &            /(tanh(beta_H*eta_H)+tanh(beta_H*(1.0-eta_H)))

        dHSVector(i)=beta_H*(1.0-(tanh(beta_H*(x-eta_H)))**2.0)
     &            /(tanh(beta_H*eta_H)+tanh(beta_H*(1.0-eta_H)))

 
      enddo

c       write(500,*) FilterMatrix,"...",colFilter,rowFilter

c      close(500)
!
      return
      end
