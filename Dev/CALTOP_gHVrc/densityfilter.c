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

void densityfilter(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
	                ITG *ne,
                  double *ttime, double *timepar,
	                ITG *mortar,double *rmin,ITG *filternnz,
                  double *FilterMatrixs,ITG *rowFilters,ITG *colFilters,ITG *filternnzElems, ITG itertop,ITG *fnnzassumed){

  char *lakon=NULL;

  ITG i,j,ne0,*kon=NULL,*ipkon=NULL;

  double *tper;
  

  double dtime,time;

  ipkon=*ipkonp;lakon=*lakonp;
  kon=*konp;

  tper=&timepar[1];

  time=*tper;
  dtime=*tper;

  ne0=*ne;


  if(itertop==1){
    
              double *elCentroid=NULL; //pointer to store Centroid of elements
              NNEW(elCentroid,double,3*ne0);  //allocate memory to element CG, initialize to 0 



              // calculate Centroid of elements
              mafillsmmain_filter(co,nk,kon,ipkon,lakon,ne,ttime,&time,mortar,&ne0,elCentroid);
              printf("\n ...Element centroids calculated... \n");

        /*      FILE *elCentroid_file;
         *     elCentroid_file=fopen("Centroidindensityfilter.dat","w"); //open in write mode
         *     for(int iii=0;iii<3*ne0;iii++){
         *               fprintf(elCentroid_file,"%.3f\n",elCentroid[iii]);
         *               }
         *     fclose(elCentroid_file);
         */

          // calculate Density Filter for elements
            mafillsmmain_filter2(ipkon,rmin,filternnz,ne,ttime,&time,&ne0,elCentroid,
                                    FilterMatrixs,rowFilters,colFilters,filternnzElems,fnnzassumed);
            printf("\n ...Distance matrix calculated... \n");
            SFREE(elCentroid);


              FILE *drow; FILE *dcol; FILE *dnnz; FILE *dval;			

/*
//////////////////////

              //Write number of non zero filter values for each element
              dnnz=fopen("dnnznotexpanded.dat","w"); //open in write mode
              for(int iii=0;iii<ne0;iii++){
                              fprintf(dnnz,"%d\n",filternnzElems[iii]);
                  }
              fclose(dnnz);




/////////
*/






            //Go through each nnz and copy to respective other half, must be in serial
            FORTRAN(mafillsm_expandfilter,(FilterMatrixs,filternnzElems,rowFilters,colFilters,ne,ttime,&time,&ne0,fnnzassumed));
            printf("\n ...Filter matrix expanded... \n");


   //           FILE *drow; FILE *dcol; FILE *dnnz; FILE *dval; 
              
              //Write non zero row values for density filter
              drow=fopen("drow.dat","w"); //open in write mode
              for(int iii=0;iii< (*fnnzassumed)*(ne0);iii++){
                        if(FilterMatrixs[iii]>0){
                              fprintf(drow,"%d\n",rowFilters[iii]);
                          }
                  }
              fclose(drow);

              //Write non zero col values for density filter
              dcol=fopen("dcol.dat","w"); //open in write mode
              for(int iii=0;iii<(*fnnzassumed)*ne0;iii++){
                        if(FilterMatrixs[iii]>0){
                              fprintf(dcol,"%d\n",colFilters[iii]);
                          }
                  }
              fclose(dcol);

              //Write non zero filter values for density filter
              dval=fopen("dval.dat","w"); //open in write mode
              for(int iii=0;iii<*fnnzassumed * ne0;iii++){
                        if(FilterMatrixs[iii]>0){
                              fprintf(dval,"%.6f\n",FilterMatrixs[iii]);
                          }
                  }
              fclose(dval);



/*
//////////////////////

              //Write number of non zero filter values for each element
              dnnz=fopen("dnnzexpanded.dat","w"); //open in write mode
              for(int iii=0;iii<ne0;iii++){
                              fprintf(dnnz,"%d\n",filternnzElems[iii]);
                  }
              fclose(dnnz);




/////////
		
*/



              //Write number of non zero filter values for each element
              dnnz=fopen("dnnz.dat","w"); //open in write mode
              for(int iii=0;iii<ne0;iii++){
                              fprintf(dnnz,"%d\n",filternnzElems[iii]);
                              if(filternnzElems[iii]>*fnnzassumed){
                                printf("\n ******* \n *************** \n *WARNING: Number of elements,%d inside filter radius for element %d exceeds %d \n *************\n **********\n ",
                                                        filternnzElems[iii],iii,*fnnzassumed);
                                printf(" \n *CAUTION: ADJUST FILTER SETTING \n");
                                exit(1);
                              }
                  }
              fclose(dnnz);
/*To write a consolidated filter file
	          FILE *filter_file2;
            filter_file2=fopen("filter.dat","w"); //open in write mode
            for(int iii=0;iii<*fnnzassumed * ne0;iii++){
                  if(FilterMatrixs[iii]>0){
  	                fprintf(filter_file2,"%d , %d , %.15f\n",rowFilters[iii],colFilters[iii],FilterMatrixs[iii]);
	                  }
              }
            fclose(filter_file2);          


            printf("\n ...Density filter written to a file... \n");
*/

      }else{

        printf("Reading density filter matrix from file");

         //Read non zeros in each row from dnnz.dat and calculate total nnzs
              *filternnz=0; //initialize
              int iii=0;
		FILE *dnnzw;
              dnnzw=fopen("dnnz.dat","r"); //open in read mode

              if (dnnzw!=NULL){
                for (iii=0;iii<ne0;iii++)
                {
                    fscanf(dnnzw,"%d",&filternnzElems[iii]);
                    *filternnz+=filternnzElems[iii];
                              if(filternnzElems[iii]>*fnnzassumed){
                                printf("\n ******* \n *************** \n *WARNING: Number of elements,%d inside filter radius for element %d exceeds %d \n *************\n **********\n ",
                                                        filternnzElems[iii],iii,*fnnzassumed);
                                printf("\n  *CAUTION: ADJUST FILTER SETTING \n");
                                exit(1);
                              }
                                    }

            fclose(dnnzw);
            }else{
                    perror("Error reading dval.dat");
                    exit(1);
            }


            double *dval=NULL; //pointer to store density filter values
            NNEW(dval,double,*filternnz);  //allocate memory to dval 

            ITG *drow=NULL; //pointer to store density filter rows 
            NNEW(drow,ITG,*filternnz);  

            ITG *dcol=NULL; //pointer to store density filter rows 
            NNEW(dcol,ITG,*filternnz);  

            FILE *dcolw;
            dcolw=fopen("dcol.dat","r"); //open in read mode

              if (dcolw!=NULL){
                for (iii=0;iii<*filternnz;iii++)
                {
                    fscanf(dcolw,"%d",&dcol[iii]);
                }

            fclose(dcolw);
            }else{
                    perror("Error reading dcol.dat");
            }

            FILE *droww;
            droww=fopen("drow.dat","r"); //open in read mode

              if (droww!=NULL){
                for (iii=0;iii<*filternnz;iii++)
                {
                    fscanf(droww,"%d",&drow[iii]);
                }

            fclose(droww);
            }else{
                    perror("Error reading drow.dat");
            }

            FILE *dvalw;
            dvalw=fopen("dval.dat","r"); //open in read mode

              if (dvalw!=NULL){
                for (iii=0;iii<*filternnz;iii++)
                {
                    fscanf(dvalw,"%lf",&dval[iii]);
                }

            fclose(dvalw);
            }else{
                    perror("Error reading dval.dat");
            }


            FORTRAN(readfilter,(FilterMatrixs,filternnzElems,rowFilters,colFilters,ne,ttime,&time,&ne0,filternnz,drow,dcol,dval,fnnzassumed));
            printf("\n...Density filter loaded from file...\n");

            SFREE(dcol);SFREE(drow);SFREE(dval);
          }
  (*ttime)+=(*tper);
  return;
}
