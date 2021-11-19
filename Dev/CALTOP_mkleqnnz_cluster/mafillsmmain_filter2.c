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


static ITG *ne1,*ne01,num_cpus,*neapar=NULL,*nebpar=NULL,*filternnz1,
            *ipkon1,
            *filternnzElems1,*rowFilters1,*colFilters1, *elarr=NULL,*fnnzassumed1;

static double *ttime1,*time1,*elCentroid1,*FilterMatrixs1,*rmin1;

void mafillsmmain_filter2(ITG *ipkon,double *rmin,ITG *filternnz,
	       ITG *ne,double *ttime,double *time,
	       ITG *ne0,double *elCentroid,double *FilterMatrixs, ITG *rowFilters,ITG *colFilters,ITG *filternnzElems,ITG *fnnzassumed ){

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
    //num_cpus=1;// overwrite for now to 1 CPU only

    if(*ne<num_cpus) num_cpus=*ne;

    pthread_t tid[num_cpus];

    /* determining the element bounds in each thread */

    NNEW(neapar,ITG,num_cpus);
    NNEW(nebpar,ITG,num_cpus);
    
    NNEW(elarr,ITG,*ne);
	//elementcpuload(neapar,nebpar,ne0,ipkon,&num_cpus);
    elementcpuload_density(neapar,nebpar,ne0,ipkon,&num_cpus,elarr);

//////////////////////////////////////////////////////////////////////////
    FILE *f1;
    f1=fopen("elemcpuarrange.dat","w"); //open in write mode
    int iii;  //counter
    fprintf(f1,"\n \n neapar \t nebpar \n");
        for (iii=0;iii<num_cpus;iii++)
            {
                fprintf(f1,"%d \t %d \n",neapar[iii],nebpar[iii]);
            }


    fclose(f1);





    ipkon1=ipkon;ne1=ne;
    ttime1=ttime;time1=time;ne01=ne0;
    elCentroid1=elCentroid;
    
    FilterMatrixs1=FilterMatrixs;
    rmin1=rmin;
    fnnzassumed1=fnnzassumed;
    filternnz1=filternnz;

  filternnzElems1=filternnzElems;
    rowFilters1=rowFilters;
    colFilters1=colFilters;

   

    printf(" Using up to %" ITGFORMAT " cpu(s) for the filter matrix calculation.\n\n", num_cpus);

    /* create threads and wait */

    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafillsmfilter2mt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

    SFREE(ithread);SFREE(neapar);SFREE(nebpar);SFREE(elarr);

    /*      for(i=0;i<num_cpus;i++){
      for(k=i*neq[1];k<i*neq[1]+neq[1];++k){printf("fext=%" ITGFORMAT ",%f\n",k-i*neq[1],fext1[k]);}
      for(k=i*neq[1];k<i*neq[1]+neq[1];++k){printf("ad=%" ITGFORMAT ",%f\n",k-i*neq[1],ad1[k]);}
      for(k=i*nzs[2];k<i*nzs[2]+nzs[2];++k){printf("au=%" ITGFORMAT ",%f\n",k-i*nzs[2],au1[k]);}
      }*/






  return;

}

/* subroutine for multithreading of mafillsm */

void *mafillsmfilter2mt(ITG *i){

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




    FORTRAN(mafillsm_filter2,(ne1,ttime1,time1,ne01,&nea,&neb,
                              elCentroid1,rmin1,
                              filternnz1,
                              FilterMatrixs1,rowFilters1,colFilters1,filternnzElems1,elarr,fnnzassumed1));

    return NULL;
}
