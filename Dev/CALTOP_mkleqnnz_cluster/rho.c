#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

void rho(double *design,int ne){
/*Reads rho values from a file density.dat
INPUT: 1)design= a vector to store rho values
       2)ne=number of elements*/
FILE *rhoFile;

    rhoFile=fopen("density.dat","r"); //open in read mode
    int i;  //counter

    //Initialize rho=1
    for(i=0;i<ne;i++){
        design[i]=1.0000;
    }

    if(rhoFile==NULL) //if rhoFile doesn't exist already,create
    {
        printf("\n...density.dat not found, initialized to 1");
        rhoFile=fopen("density.dat","w");
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
        }

        fclose(rhoFile);
}
