#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

void rho(double *design,int ne, double pstiff){
/*Reads rho values from a file densityfsi.dat
INPUT: 1)design= a vector to store rho values
       2)ne=number of elements
       3) pstiff flag*/
FILE *rhoFile;


    int i;  //counter

    //Initialize rho=1
    for(i=0;i<ne;i++){
        design[i]=1.0000;
    }

    if (pstiff!=0){  
        rhoFile=fopen("densityfsi.dat","r"); //open in read mode
        if(rhoFile==NULL) //if rhoFile doesn't exist already,create
        {
            printf("\nNOTE: densityfsi.dat not found, initialized to 1 . \n");
            rhoFile=fopen("densityfsi.dat","w");
            for (i=0;i<ne;i++)
                {
                    fprintf(rhoFile,"%.15f \n",design[i]);
                }
        }

        else{//if file exists,           
            for (i=0;i<ne;i++)
                {
                    fscanf(rhoFile,"%lf",&design[i]);
                }
            printf("\nFiltered rho values read from densityfsi.dat .\n"); 
            }

        fclose(rhoFile);
    }   //end: if pstiff !=0
}
