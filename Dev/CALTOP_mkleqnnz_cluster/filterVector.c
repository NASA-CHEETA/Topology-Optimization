/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"

void filterVector(ITG **ipkonp,double *Vector, double *VectorFiltered,double *FilterMatrix,
ITG *filternnzElem,ITG *rowFilter,
ITG *colFilter,ITG *ne,double *ttime, double *timepar, ITG *fnnzassumed, double *q){



  ITG i,j,ne0,*ipkon=NULL;

  


  double *tper;
  
 

  double dtime,time;

  ipkon=*ipkonp;
  tper=&timepar[1];
  time=*tper;
  dtime=*tper;
  ne0=*ne;
   

/*      FILE *elCentroid_file;
       elCentroid_file=fopen("Centroidindensityfilter.dat","w"); //open in write mode
      for(int iii=0;iii<3*ne0;iii++){
        fprintf(elCentroid_file,"%.3f\n",elCentroid[iii]);
        }
        fclose(elCentroid_file);
*/


// calculate Apply filter to a vector
  mafillsmmain_Vectorfilter(ipkon,Vector,VectorFiltered,FilterMatrix,filternnzElem,
            rowFilter,colFilter,ne,ttime,&time,&ne0, fnnzassumed,q);


/*      FILE *filter_file;
       filter_file=fopen("Filterfilte.dat","w"); //open in write mode
      for(int iii=0;iii<100*ne0;iii++){
        fprintf(filter_file,"%d , %d , %.3f\n",rowFilter[iii],colFilter[iii],FilterMatrix[iii]);
        }
        fclose(filter_file);

*/




  (*ttime)+=(*tper);

  return;
}
