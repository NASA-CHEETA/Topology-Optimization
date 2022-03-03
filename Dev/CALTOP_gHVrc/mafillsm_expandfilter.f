      subroutine mafillsm_expandfilter(FilterMatrixs,filternnzElems,
     &  rowFilters,colFilters,
     &  ne,ttime,time,ne0,fnnzassumed)
!
!     expand the symmetric FilterMatrixs
!
      implicit none
!
c      integer 
!
!
      integer filternnzElems(*),ne0,filternnz,
     & fnnzassumed,
     &  ne,i,j,rowFilters(fnnzassumed,*),colFilters(fnnzassumed,*),
     &  rowval,colval
!
      real*8 FilterMatrixs(fnnzassumed,*),ttime,
     &  time,value
!
      intent(in) ne,
     &  ttime,time,
     &  ne0, fnnzassumed
!
      intent(inout) FilterMatrixs,
     &  filternnzElems,
     &  rowFilters,colFilters

!
c      open (unit=200,file="filternnzelem.dat")
 
!s
!       
        
!      loop over all elements
!
      do i=1,ne
        do j=1,filternnzElems(i)
                  rowval=rowFilters(j,i)
                  colval=colFilters(j,i)
                  value=FilterMatrixs(j,i)


                

          if(colval.GT.i) then
            
            filternnzElems(colval)=filternnzElems(colval)+1

            rowFilters(filternnzElems(colval),colval)=colval
            colFilters(filternnzElems(colval),colval)=i
            FilterMatrixs(filternnzElems(colval),colval)=value
          endif

        enddo

c	write(200,*) i,filternnzElems(i)
      enddo

c       write(200,*) FilterMatrix,"...",colFilter,rowFilter

c        close(200)
!

      return
      end
