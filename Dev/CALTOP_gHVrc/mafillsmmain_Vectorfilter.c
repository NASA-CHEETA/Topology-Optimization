/*     CalculiX - A 3-dimensional finite element program                 */
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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"


static ITG *ne1,*ne01,num_cpus,*neapar=NULL,*nebpar=NULL,
            *ipkon1,*filternnzElem1,*rowFilter1,*colFilter1, *fnnzassumed1;

static double *ttime1,*time1,*FilterMatrix1,*Vector1,*VectorFiltered1, *q1;

void mafillsmmain_Vectorfilter(ITG *ipkon,double *Vector,double *VectorFiltered,
    double *FilterMatrix,ITG *filternnzElem,
           ITG *rowFilter,ITG *colFilter,
	       ITG *ne,double *ttime,double *time,
	       ITG *ne0, ITG *fnnzassumed, double *q){

    ITG i,j;

    /* variables for multithreading procedure */

    ITG sys_cpus,*ithread=NULL;
    char *env,*envloc,*envsys;

    num_cpus = 0;
    sys_cpus=0;

    /* explicit user declaration prevails */

    envsys=getenv("NUMBER_OF_CPUS");
    if(envsys){
	sys_cpus=atoi(envsys);
	if(sys_cpus<0) sys_cpus=0;
    }

//    sys_cpus=1;

    /* automatic detection of available number of processors */

    if(sys_cpus==0){
	sys_cpus = getSystemCPUs();
	if(sys_cpus<1) sys_cpus=1;
    }

    /* local declaration prevails, if strictly positive */

    envloc = getenv("CCX_NPROC_STIFFNESS");
    if(envloc){
	num_cpus=atoi(envloc);
	if(num_cpus<0){
	    num_cpus=0;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}

    }

    /* else global declaration, if any, applies */

    env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
	if (env)
	    num_cpus = atoi(env);
	if (num_cpus < 1) {
	    num_cpus=1;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
    }

// next line is to be inserted in a similar way for all other paralell parts

// num_cpus=1;// overwrite for now to 1 CPU only


    if(*ne<num_cpus) num_cpus=*ne;

    pthread_t tid[num_cpus];

    /* determining the element bounds in each thread */

    NNEW(neapar,ITG,num_cpus);
    NNEW(nebpar,ITG,num_cpus);
    

	elementcpuload(neapar,nebpar,ne0,ipkon,&num_cpus);


    ipkon1=ipkon;ne1=ne; q1=q;
    ttime1=ttime;time1=time;ne01=ne0;
    
    FilterMatrix1=FilterMatrix;
    Vector1=Vector;
    VectorFiltered1=VectorFiltered;
    filternnzElem1=filternnzElem;
    rowFilter1=rowFilter;
    colFilter1=colFilter;
    fnnzassumed1=fnnzassumed;

    /* FilterVector */

    printf(" Using up to %" ITGFORMAT " cpu(s) for the filtering of vector.\n\n", num_cpus);

    /* create threads and wait */

    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafillsmVectorfiltermt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

    SFREE(ithread);SFREE(neapar);SFREE(nebpar);

    /*      for(i=0;i<num_cpus;i++){
      for(k=i*neq[1];k<i*neq[1]+neq[1];++k){printf("fext=%" ITGFORMAT ",%f\n",k-i*neq[1],fext1[k]);}
      for(k=i*neq[1];k<i*neq[1]+neq[1];++k){printf("ad=%" ITGFORMAT ",%f\n",k-i*neq[1],ad1[k]);}
      for(k=i*nzs[2];k<i*nzs[2]+nzs[2];++k){printf("au=%" ITGFORMAT ",%f\n",k-i*nzs[2],au1[k]);}
      }*/






  return;

}

/* subroutine for multithreading of mafillsm */

void *mafillsmVectorfiltermt(ITG *i){

    ITG nea,neb;

    nea=neapar[*i]+1;
    neb=nebpar[*i]+1;

/*FILE *rhoFile;

    rhoFile=fopen("densityout_mafillsmmain.dat","w"); //open in write mode


    int i;  //counter
        for (i=0;i<ne;i++)
            {
                fprintf(rhoFile,"%f \n",design1[i]);
            }
        fclose(rhoFile);*/




    FORTRAN(mafillsmvectorfilter,(ne1,ttime1,time1,ne01,&nea,&neb,
                              FilterMatrix1,Vector1,VectorFiltered1,
                              filternnzElem1,rowFilter1,colFilter1,fnnzassumed1,q1));

    return NULL;
}
