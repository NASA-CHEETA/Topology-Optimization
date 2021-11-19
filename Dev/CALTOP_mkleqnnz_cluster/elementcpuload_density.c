/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

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
#include <string.h>

#include "CalculiX.h"

void elementcpuload_density(ITG *neapar,ITG *nebpar,ITG *ne,ITG *ipkon,ITG *num_cpus,ITG *elarr){

 /*  divides the elements into ranges with an equal number of
     active elements (element numbering may have gaps) for
     parallel processing on different cpus */

    ITG i,nepar,*ipar=NULL,idelta,ideltac,ideltaf,isum;
    
    NNEW(ipar,ITG,*ne);

    nepar=0;
    for(i=0;i<*ne;i++){
	if(ipkon[i]>-1){

	    /* active element */
	    
	    ipar[nepar]=i;
	    nepar++;
	}
    }

    if(nepar<*num_cpus) *num_cpus=nepar;

    /* dividing the element number range into num_cpus equal numbers of 
       active elements */



    ideltaf=nepar/(*num_cpus); //division floor 
    ideltac=ideltaf+1;  //division ceil
        //determine whether to use ceil or floor
        if(ideltac*(*num_cpus-1) < nepar || ideltac*(*num_cpus-1) == nepar){
            idelta=ideltac;
                }else{
                idelta=ideltaf;
        }

	//If perfectly divisible, use floor always
	if(nepar%(*num_cpus)==0){
		idelta=ideltaf;
		}

    isum=0;
    for(i=0;i<*num_cpus;i++){
	neapar[i]=ipar[isum];
	if(i!=*num_cpus-1){
	    isum+=idelta;
	}else{
	    isum=nepar;
	}
	nebpar[i]=ipar[isum-1];
    }
    
        //arrange job 1st step
    int index=0;
    for(i=0;i<nepar;i+=2){
       // printf("%d",a[index]);
        elarr[i]=ipar[index];
        index++;
        
    }

        //arrange job second step
    index=nepar-1;
    for(i=1;i<nepar;i+=2){
        elarr[i]=ipar[index];
        index--;
    }
/*
     FILE *f1;

    f1=fopen("elemcpuarranged.dat","w"); //open in write mode


    int iii;  //counter
    fprintf(f1,"elarr contents\n");
        for (iii=0;iii<nepar;iii++)
            {
                fprintf(f1,"%d \n",elarr[iii]);
            }
    fprintf(f1,"\n \n neapar \t nebpar \n");
        for (iii=0;iii<*num_cpus;iii++)
            {
                fprintf(f1,"%d \t %d \n",neapar[iii],nebpar[iii]);
            }


        fclose(f1);*/

    SFREE(ipar);
    
    return;
}
