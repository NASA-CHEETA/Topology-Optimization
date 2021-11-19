      subroutine readfilter(FilterMatrixs,filternnzElems,
     &  rowFilters,colFilters,
     &  ne,ttime,time,ne0,filternnz,drow,dcol,dval,fnnzassumed)
!
!     arrange filter vectors in array
!
      implicit none
!
!
!
      integer filternnzElems(*),ne0,filternnz,fnnzassumed,
     &  ne,i,j,rowFilters(fnnzassumed,*),
     &  colFilters(fnnzassumed,*),
     &  rowval,colval,drow(*),dcol(*), index
!
      real*8 FilterMatrixs(fnnzassumed,*),ttime,
     &  time,value, dval(*)
!
      intent(in) ne,filternnz,drow,dcol,dval,
     &  ttime,time, fnnzassumed,
     &  ne0
!
      intent(inout) FilterMatrixs,
     &  rowFilters,colFilters

!
c       open (unit=300,file="readfilterf.dat")
 
!s
!       
       index=1 
!      loop over all nonzeros
!
       do i=1,filternnz
        rowval=drow(i)
        colval=dcol(i)
        value=dval(i)

        rowFilters(index,rowval)=rowval
        colFilters(index,rowval)=colval
        FilterMatrixs(index,rowval)=value

c        write(300,*) rowFilters(index,rowval),
c     &              colFilters(index,rowval),
c     &              FilterMatrixs(index,rowval)
              
        if(index.LT.filternnzElems(rowval)) then
            
            index=index + 1

        else
            index = 1 
        
        endif

       enddo


c       write(200,*) FilterMatrix,"...",colFilter,rowFilter

c        close(300)
!

      return
      end
