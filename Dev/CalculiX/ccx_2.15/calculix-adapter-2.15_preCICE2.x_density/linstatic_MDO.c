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
#include "omp.h"
#include "CalculiX.h"
#include <time.h>

/*---Adapter header declaration---*/
#include "adapter/PreciceInterface.h"

#ifdef SPOOLES
   #include "spooles.h"
#endif
#ifdef SGI
   #include "sgi.h"
#endif
#ifdef TAUCS
   #include "tau.h"
#endif
#ifdef PARDISO
   #include "pardiso.h"
#endif



/**********************************************************************/
/**
 * @brief gkdas2: Get the ith set name from SET list
 * 
 * @param string 
 * @param position 
 * @param length 
 * @return *char 
 */
char *getSubstring2(char *string, int position, int length)
{
   char *p;
   int c;
 
   p = malloc(length+1);
   if (p == NULL)
   {
      printf("Unable to allocate memory in getSubstring.\n");
      exit(1);
   }
 
   for (c = 0; c < length; c++)
   {
      *(p+c) = *(string+position-1);      
      string++;  
   }
   *(p+c) = '\0';
   return p;
}

/**
 * @brief gkdas2: write surface node set NSurface displacements
 * 
 * @param nset 
 * @param ialset 
 * @param set 
 * @param istartset 
 * @param iendset 
 */
void writeSurfNodes(ITG nset, ITG *ialset, char *set, ITG *istartset, 
					ITG *iendset, double *v)
{
	char surfaceSet[] = "NSURFACEN"; // ccx converts all strings to uppercase, adds N for nodes
	char *dummyNameset;
	ITG istart;
	ITG iend;
	ITG nodeNumber;
	
	for(int i = 0; i< nset; i++)
	{	
		dummyNameset = getSubstring2(set, i*81+1,strlen(surfaceSet));
		if (strcmp(dummyNameset, surfaceSet) == 0)
		{	
			printf("For node set: %s\n", dummyNameset);
			istart = istartset[i]-1;
			iend = iendset[i]-1;

			for (ITG iset = istart; iset <= iend;iset++)
			{	
				nodeNumber = ialset[iset];
				printf("Disp: Node: %d, UX: %.9lf, UY: %.9lf, UZ: %.9lf \n", nodeNumber, 
							v[4*(nodeNumber-1)+1], v[4*(nodeNumber-1)+2], v[4*(nodeNumber-1)+3]); // nodeNumber is 1 based
			}			
		}	
	}
	SFREE(dummyNameset);
}

/**
 * @brief gkdas2 Add displacement v to updated coordinate list. 
 * TO DO: Verify for 2D
 * @param coUpdated 
 * @param v  : latest displacement vector
 * @param nk : number of nodes
 * @param mt : dof index. 1 temperature + 3 displacement
 */
void updateCO(double *coUpdated, double *v, int nk, int mt)
{	
	int i, a, b;
	int num_cpus=0;
	char *env;

	// get num of threads declaration, if any
	env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
    	if (env)
      		num_cpus = atoi(env);
    	if (num_cpus < 1) {
      		num_cpus=1;
    	}
    }

	#pragma omp parallel for num_threads(num_cpus) private(a,b)
	for(i = 0; i<= nk-1; i++)
	{
		a = 3*i;
		b = a + i + 1;
		coUpdated[a]+=v[b];
		coUpdated[a+1]+=v[b+1];
		coUpdated[a+2]+=v[b+2];
	}
}

/**
 * @brief gkdas2  Write the updated coordinate list to updatedCOORD.csv. O
 * TO DO: Verify for 2D
 * @param coUpdated 
 * @param nk 
 * @param mt 
 */
void writeUpdatedCO(double *coUpdated, int nk, int mt)
{
    FILE *fp;

    fp = fopen("updatedCOORD.csv", "w+");
	
	fprintf(fp,"NodeID,x,y,z \n");
    for(int i = 0; i <= nk-1; i++)
    {
        fprintf(fp, "%d ,%0.6lf, %0.6lf, %0.6lf \n", i+1, coUpdated[3*i], coUpdated[3*i+1], coUpdated[3*i+2]);
    }
    fclose(fp);
}

void writeBaselineCO(double *co, int nk, int mt)
{
    FILE *fp;

    fp = fopen("BaselineCOORD.csv", "w+");
	
	fprintf(fp,"NodeID,x,y,z \n");
    for(int i = 0; i <= nk-1; i++)
    {
        fprintf(fp, "%d ,%0.6lf, %0.6lf, %0.6lf \n", i+1, co[3*i], co[3*i+1], co[3*i+2]);
    }
    fclose(fp);
}

/* gkdas2: linstaic_MDO begins here  */
void linstatic_MDO(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
	     ITG *ne, 
	     ITG *nodeboun, ITG *ndirboun, double *xboun, ITG *nboun, 
	     ITG *ipompc, ITG *nodempc, double *coefmpc, char *labmpc,
             ITG *nmpc, 
	     ITG *nodeforc, ITG *ndirforc,double *xforc, ITG *nforc, 
	     ITG *nelemload, char *sideload, double *xload,
	     ITG *nload, ITG *nactdof, 
	     ITG **icolp, ITG *jq, ITG **irowp, ITG *neq, ITG *nzl, 
	     ITG *nmethod, ITG *ikmpc, ITG *ilmpc, ITG *ikboun, 
	     ITG *ilboun,
	     double *elcon, ITG *nelcon, double *rhcon, ITG *nrhcon,
	     double *alcon, ITG *nalcon, double *alzero, ITG **ielmatp,
	     ITG **ielorienp, ITG *norien, double *orab, ITG *ntmat_,
	     double *t0, double *t1, double *t1old,
	     ITG *ithermal,double *prestr, ITG *iprestr, 
	     double *vold,ITG *iperturb, double *sti, ITG *nzs,  
	     ITG *kode, char *filab, double *eme,
             ITG *iexpl, double *plicon, ITG *nplicon, double *plkcon,
             ITG *nplkcon,
             double **xstatep, ITG *npmat_, char *matname, ITG *isolver,
             ITG *mi, ITG *ncmat_, ITG *nstate_, double *cs, ITG *mcs,
             ITG *nkon, double **enerp, double *xbounold,
	     double *xforcold, double *xloadold,
             char *amname, double *amta, ITG *namta,
	     ITG *nam, ITG *iamforc, ITG *iamload,
             ITG *iamt1, ITG *iamboun, double *ttime, char *output, 
             char *set, ITG *nset, ITG *istartset,
             ITG *iendset, ITG *ialset, ITG *nprint, char *prlab,
             char *prset, ITG *nener, double *trab, 
             ITG *inotr, ITG *ntrans, double *fmpc, char *cbody, ITG *ibody,
	     double *xbody, ITG *nbody, double *xbodyold, double *timepar,
	     double *thicke, char *jobnamec,char *tieset,ITG *ntie,
	     ITG *istep,ITG *nmat,ITG *ielprop,double *prop,char *typeboun,
	     ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,ITG *icontact,
             char *orname,double *design, double *penal, char *preciceParticipantName, char *configFilename,ITG *ikforc,ITG *ilforc)
{
  
  char description[13]="            ",*lakon=NULL,stiffmatrix[132]="",
       fneig[132]="",jobnamef[396]="";

  ITG *inum=NULL,k,*icol=NULL,*irow=NULL,ielas=0,icmd=0,iinc=1,nasym=0,i,j,ic,ir,
      mass[2]={0,0}, stiffness=1, buckling=0, rhsi=1, intscheme=0,*ncocon=NULL,
      *nshcon=NULL,mode=-1,noddiam=-1,*ipobody=NULL,inewton=0,coriolis=0,iout,
      ifreebody,*itg=NULL,ntg=0,symmetryflag=0,inputformat=0,ngraph=1,im,
      mt=mi[1]+1,ne0,*integerglob=NULL,iglob=0,*ipneigh=NULL,*neigh=NULL,
      icfd=0,*inomat=NULL,*islavact=NULL,*islavnode=NULL,*nslavnode=NULL,
      *islavsurf=NULL,nretain,*iretain=NULL,*noderetain=NULL,*ndirretain=NULL,
      nmethodl,nintpoint,ifacecount,memmpc_,mpcfree,icascade,maxlenmpc,
      ncont=0,*itietri=NULL,*koncont=NULL,nslavs=0,ismallsliding=0,
      *itiefac=NULL,*imastnode=NULL,*nmastnode=NULL,*imastop=NULL,iitsta,
      *iponoels=NULL,*inoels=NULL,*ipe=NULL,*ime=NULL,iit=-1,iflagact=0,
      icutb=0,*kon=NULL,*ipkon=NULL,*ielmat=NULL,ialeatoric=0,kscale=1,
      *iponoel=NULL,*inoel=NULL,zero=0,nherm=1,nev=*nforc,node,idir,
      *ielorien=NULL,network=0,nrhs=1,iperturbsav;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*xstiff=NULL,*stiini=NULL,*tper,
         *f=NULL,*fn=NULL,qa[4],*fext=NULL,*epn=NULL,*xstateini=NULL,
         *vini=NULL,*stx=NULL,*enern=NULL,*xbounact=NULL,*xforcact=NULL,
         *xloadact=NULL,*t1act=NULL,*ampli=NULL,*xstaten=NULL,*eei=NULL,
         *enerini=NULL,*cocon=NULL,*shcon=NULL,*physcon=NULL,*qfx=NULL,
         *qfn=NULL,sigma=0.,*cgr=NULL,*xbodyact=NULL,*vr=NULL,*vi=NULL,
         *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*springarea=NULL,
         *eenmax=NULL,*fnr=NULL,*fni=NULL,*emn=NULL,*clearini=NULL,ptime,
         *emeini=NULL,*doubleglob=NULL,*au=NULL,*ad=NULL,*b=NULL,*aub=NULL,
         *adb=NULL,*pslavsurf=NULL,*pmastsurf=NULL,*cdn=NULL,*cdnr=NULL,
         *cdni=NULL,*submatrix=NULL,*xnoels=NULL,*cg=NULL,*straight=NULL,
         *areaslav=NULL,*xmastnor=NULL,theta=0.,*ener=NULL,*xstate=NULL,
         *fnext=NULL,*energyini=NULL,*energy=NULL,*d=NULL,alea=0.1;

  FILE *f1,*f2;
  
 #ifdef SGI
  ITG token;
 #endif

	/** @gkdas2
	 * Create a new vector to store updated coordinates and copy initial coordinates
	 */
	double *coUpdated = NULL;
	NNEW(coUpdated,double,3**nk);
	memcpy(coUpdated, co, 3**nk*sizeof(double));




  /* dummy arguments for the results call */
  double *veold=NULL,*accold=NULL,bet,gam,dtime,time,reltime=1.;

  irow=*irowp;
  ener=*enerp;
  xstate=*xstatep;
  ipkon=*ipkonp;
  lakon=*lakonp;
  kon=*konp;
  ielmat=*ielmatp;
  ielorien=*ielorienp;
  icol=*icolp;
  
  for(k=0;k<3;k++)
  {
    strcpy1(&jobnamef[k*132],&jobnamec[k*132],132);
  }

  tper=&timepar[1];

  time=*tper;
  dtime=*tper;

  ne0=*ne;

  /*---Asign values to CCX structure to hold coupling variables---*/
  struct SimulationData simulationData = 
  {
    .ialset = ialset,
    .ielmat = ielmat,
    .istartset = istartset,
    .iendset = iendset,
    .kon = kon,
    .ipkon = ipkon,
    .lakon = &lakon,
    .co = co,
    .set = set,
    .nset = *nset,
    .ikboun = ikboun,
    .ikforc = ikforc,
    .ilboun = ilboun,
    .ilforc = ilforc,
    .nboun = *nboun,
    .nforc = *nforc,
    .nelemload = nelemload,
    .nload = *nload,
    .sideload = sideload,
    .mt = mt,
    .nk = *nk,
    .theta = &theta,
    //.dtheta = &dtheta,
    .tper = tper,
    .nmethod = nmethod,
    .xload = xload,
    .xforc = xforc,
    .xboun = xboun,
    .ntmat_ = ntmat_,
    .vold = vold,
    .veold = veold,
    .fn = fn,
    .cocon = cocon,
    .ncocon = ncocon,
    .mi = mi
  };





	/*---determining the global values to be used as boundary conditions
	for a submodel 

	iglob=-1 if global results are from a *FREQUENCY calculation
	iglob=0 if no global results are used by boundary conditions
	iglob=1 if global results are from a *STATIC calculation ---*/ 

	getglobalresults(jobnamec,&integerglob,&doubleglob,nboun,iamboun,xboun,
		   nload,sideload,iamload,&iglob,nforc,iamforc,xforc,
                   ithermal,nk,t1,iamt1,&sigma);

  	/*---allocating fields for the actual external loading---*/
	NNEW(xbounact,double,*nboun);

  	for(k=0;k<*nboun;++k)
  	{
		xbounact[k]=xbounold[k];
  	}

  	NNEW(xforcact,double,*nforc);

  	NNEW(xloadact,double,2**nload);

  	NNEW(xbodyact,double,7**nbody);

  	/*---copying the rotation axis and/or acceleration vector---*/ 
  	for(k=0;k<7**nbody;k++)
  	{
		xbodyact[k]=xbody[k];
  	}

  	/*---assigning the body forces to the elements---*/ 
  	if(*nbody>0)
  	{
    	ifreebody=*ne+1;

    	NNEW(ipobody,ITG,2*ifreebody**nbody);

    	for(k=1;k<=*nbody;k++)
		{
			FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			     iendset,ialset,&inewton,nset,&ifreebody,&k));

	  		RENEW(ipobody,ITG,2*(*ne+ifreebody));
    	}

      RENEW(ipobody,ITG,2*(ifreebody-1));
  	}

    /*---Adapter: Create the interfaces and initialize the coupling---*/
	printf("About to enter Aeroelastic domain in Calculix with names %s and %s \n", preciceParticipantName, configFilename);

	Precice_Setup( configFilename, preciceParticipantName, &simulationData );

	while (Precice_IsCouplingOngoing() )
	{
	  /*---Adapter: Adjust solver time step---*/
    	Precice_AdjustSolverTimestep( &simulationData );

      /*---See if coupling data (aerodynamic tractions)need to be read---*/
	  	Precice_ReadCouplingData(&simulationData);

  	  /*---allocating a field for the instantaneous amplitude---*/
  		NNEW(ampli,double,*nam);

  		FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
	      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
          xbounold,xboun,xbounact,iamboun,nboun,
	      nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
          ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
          iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
          ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
          shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon));

  		/* determining the internal forces and the stiffness coefficients */
  		NNEW(f,double,*neq);

  		/* allocating a field for the stiffness matrix */
  		NNEW(xstiff,double,(long long)27*mi[0]**ne);

  		/* for a *STATIC,PERTURBATION analysis with submodel boundary
        conditions from a *FREQUENCY analysis iperturb[0]=1 has to be
        temporarily set to iperturb[0]=0 in order for f to be calculated in
        resultsini and subsequent results* routines */
  		if((*nmethod==1)&&(iglob<0)&&(iperturb[0]>0))
		{
      		iperturbsav=iperturb[0];
      		iperturb[0]=0;
  		}

  		iout=-1;
  		NNEW(v,double,mt**nk);
  		NNEW(fn,double,mt**nk);
  		NNEW(stx,double,6*mi[0]**ne);
  		NNEW(inum,ITG,*nk);

  		results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	  	elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	  	ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	  	prestr,iprestr,filab,eme,emn,een,iperturb,
	  	f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	  	ndirboun,xbounact,nboun,ipompc,
	  	nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,
	  	&bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	  	xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	  	&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	  	emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	  	iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	  	fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	  	&reltime,&ne0,thicke,shcon,nshcon,
	  	sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	  	mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	  	islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
          inoel,nener,orname,&network,ipobody,xbodyact,ibody,typeboun);

  		SFREE(v);SFREE(fn);SFREE(stx);SFREE(inum);
  		iout=1;

  		if((*nmethod==1)&&(iglob<0)&&(iperturb[0]>0))
		{
      		iperturb[0]=iperturbsav;
  		}
  
  		/* determining the system matrix and the external forces */

  		NNEW(ad,double,*neq);
  		NNEW(fext,double,*neq);

		/* linear static calculation */

    	NNEW(au,double,*nzs);
    	nmethodl=*nmethod;

    	/* if submodel calculation with a global model obtained by
        a *FREQUENCY calculation: replace stiffness matrix K by
        K-sigma*M */

    	if(iglob<0)
		{
	  		mass[0]=1;
	  		NNEW(adb,double,*neq);
	  		NNEW(aub,double,nzs[1]);
    	}

  		mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounact,nboun,
	    	ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	    	nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
	    	nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,&nmethodl,
	    	ikmpc,ilmpc,ikboun,ilboun,
	    	elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    	ielorien,norien,orab,ntmat_,
	    	t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
	    	nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	    	xstiff,npmat_,&dtime,matname,mi,
            ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,physcon,
            shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,&coriolis,
	    	ibody,xloadold,&reltime,veold,springarea,nstate_,
            xstateini,xstate,thicke,integerglob,doubleglob,
	    	tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
	    	pmastsurf,mortar,clearini,ielprop,prop,&ne0,fnext,&kscale,
	    	iponoel,inoel,&network,ntrans,inotr,trab,design,penal);

  		/* check for negative Jacobians */

  		if(nmethodl==0) *nmethod=0;

  		if(nasym==1)
		{
      		RENEW(au,double,2*nzs[1]);
      		symmetryflag=2;
      		inputformat=1;
      
      		FORTRAN(mafillsmas,(co,nk,kon,ipkon,lakon,ne,nodeboun,
                  ndirboun,xbounact,nboun,
		  	ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		  	nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		  	nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		  	nmethod,ikmpc,ilmpc,ikboun,ilboun,
		  	elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		  	ielmat,ielorien,norien,orab,ntmat_,
		  	t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
		  	nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		  	xstiff,npmat_,&dtime,matname,mi,
                  ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
                  physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
                  &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
                  xstateini,xstate,thicke,
                  integerglob,doubleglob,tieset,istartset,iendset,
		  	ialset,ntie,&nasym,pslavsurf,pmastsurf,mortar,clearini,
		  	ielprop,prop,&ne0,&kscale,iponoel,inoel,&network));
  		}

  		/* determining the right hand side */
  		NNEW(b,double,*neq);

  		for(k=0;k<*neq;++k)
		{
      		b[k]=fext[k]-f[k];
  		}

  		SFREE(fext);
		SFREE(f);
  		
		/*-----------------------------------------------------------------*/  
		/*-------------------------------LINEAR STATIC---------------------*/
		/*-----------------------------------------------------------------*/
  		if(*nmethod!=0)
  		{

    		/* linear static applications */

    		if(*isolver==0)
			{
				#ifdef SPOOLES
      			spooles(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,&symmetryflag,
              	&inputformat,&nzs[2]);
				#else
            	printf("*ERROR in linstatic: the SPOOLES library is not linked\n\n");
            	FORTRAN(stop,());
				#endif
    		
			}
    		else if((*isolver==2)||(*isolver==3))
			{
      			if(nasym>0)
				{
	  				printf(" *ERROR in nonlingeo: the iterative solver cannot be used for asymmetric matrices\n\n");
	  				FORTRAN(stop,());
      			}
      			preiter(ad,&au,b,&icol,&irow,neq,nzs,isolver,iperturb);
    		}
    		else if(*isolver==4)
			{
				#ifdef SGI
      			if(nasym>0)
				{
	  				printf(" *ERROR in nonlingeo: the SGI solver cannot be used for asymmetric matrices\n\n");
	  				FORTRAN(stop,());
      			}
      			token=1;
      			sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,token);
				#else
            	printf("*ERROR in linstatic: the SGI library is not linked\n\n");
            	FORTRAN(stop,());
				#endif
    		}
    		else if(*isolver==5)
			{
				#ifdef TAUCS
      			if(nasym>0)
				{
	  				printf(" *ERROR in nonlingeo: the TAUCS solver cannot be used for asymmetric matrices\n\n");
	  				FORTRAN(stop,());
      			}
      			tau(ad,&au,adb,aub,&sigma,b,icol,&irow,neq,nzs);
				#else
            	printf("*ERROR in linstatic: the TAUCS library is not linked\n\n");
            	FORTRAN(stop,());
				#endif
    		}
    		else if(*isolver==7)
			{
				#ifdef PARDISO
      			pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,
		   		&symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
				#else
            	printf("*ERROR in linstatic: the PARDISO library is not linked\n\n");
            	FORTRAN(stop,());
				#endif
    		}

    		/* saving of ad and au for sensitivity analysis */

    		for(i=0;i<*ntie;i++)
			{
				if(strcmp1(&tieset[i*243+80],"D")==0)
				{
	    
	    			strcpy(stiffmatrix,jobnamec);
	    			strcat(stiffmatrix,".stm");
	    
	    			if((f1=fopen(stiffmatrix,"wb"))==NULL)
					{
						printf("*ERROR in linstatic: cannot open stiffness matrix file for writing...");
						exit(0);
	    			}
	    
	    			/* storing the stiffness matrix */

            		/* nzs,irow,jq and icol have to be stored too, since the static analysis
               		can involve contact, whereas in the sensitivity analysis contact is not
               		taken into account while determining the structure of the stiffness
               		matrix (in mastruct.c)
	     			*/
	    
	    			if(fwrite(&nasym,sizeof(ITG),1,f1)!=1)
					{
						printf("*ERROR saving the symmetry flag to the stiffness matrix file...");
						exit(0);
	    			}
	    			if(fwrite(nzs,sizeof(ITG),3,f1)!=3)
					{
						printf("*ERROR saving the number of subdiagonal nonzeros to the stiffness matrix file...");
						exit(0);
	    			}
	    			if(fwrite(irow,sizeof(ITG),nzs[2],f1)!=nzs[2])
					{
						printf("*ERROR saving irow to the stiffness matrix file...");
						exit(0);
	    			}
	    			if(fwrite(jq,sizeof(ITG),neq[1]+1,f1)!=neq[1]+1)
					{
						printf("*ERROR saving jq to the stiffness matrix file...");
						exit(0);
	    			}
	    			if(fwrite(icol,sizeof(ITG),neq[1],f1)!=neq[1])
					{
						printf("*ERROR saving icol to the stiffness matrix file...");
						exit(0);
	    			}
	    			if(fwrite(ad,sizeof(double),neq[1],f1)!=neq[1])
					{
						printf("*ERROR saving the diagonal of the stiffness matrix to the stiffness matrix file...");
						exit(0);
	    			}

	    			if(fwrite(au,sizeof(double),nzs[2],f1)!=nzs[2])
					{
						printf("*ERROR saving the off-diagonal terms of the stiffness matrix to the stiffness matrix file...");
						exit(0);
	    			}
	    			fclose(f1);

	    			break;
        		}
    		}
    
    		SFREE(ad);
			SFREE(au);

    		if(iglob<0)
			{
				SFREE(adb);
				SFREE(aub);
			}

    		/* calculating the displacements and the stresses and storing */
    		/* the results in frd format for each valid eigenmode */

    		NNEW(v,double,mt**nk);
    		NNEW(fn,double,mt**nk);
    		NNEW(stn,double,6**nk);
    		NNEW(inum,ITG,*nk);
    		NNEW(stx,double,6*mi[0]**ne);
  
    		if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
    		if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
    		if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
    		if(strcmp1(&filab[2175],"CONT")==0) NNEW(cdn,double,6**nk);

    		NNEW(eei,double,6*mi[0]**ne);
    		if(*nener==1)
			{
				NNEW(stiini,double,6*mi[0]**ne);
        		NNEW(emeini,double,6*mi[0]**ne);
				NNEW(enerini,double,mi[0]**ne);
			}

    		results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	    	elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    	ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	    	prestr,iprestr,filab,eme,emn,een,iperturb,
            f,fn,nactdof,&iout,qa,vold,b,nodeboun,ndirboun,xbounact,nboun,ipompc,
	    	nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,&bet,
            &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    	xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
            ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    	nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
            &ne0,thicke,shcon,nshcon,
            sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
            mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	    	islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
            inoel,nener,orname,&network,ipobody,xbodyact,ibody,typeboun);

			simulationData.fn = fn;

			/*---Write the displacement data---*/


    		SFREE(eei);
    		if(*nener==1)
			{
				SFREE(stiini);
				SFREE(emeini);
				SFREE(enerini);
			}

			/*---Copy current v into vold---*/

    		memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);

			/* gkdas2: add displacement to coordinates*/
			printf("Updating coordinates to heap\n");
			updateCO(coUpdated, v, *nk, mt);

			printf("Displacement of surface nodes:\n");
			writeSurfNodes(*nset, ialset, set, istartset, iendset, v);

			/*---Save the current displacement state for implicit calculations---*/

			if(Precice_IsWriteCheckpointRequired())
    		{
    		//	Precice_WriteIterationCheckpoint( &simulationData, vold );
        		Precice_FulfilledWriteCheckpoint();
    		}

			/*---Transmit the interface discplacement state---*/ 
			Precice_WriteCouplingData(&simulationData);

			/*---Advance the coupling state---*/
			Precice_Advance(&simulationData);

			/*---Write these deformed coordinates to coUpdated---*/

	



			printf("\n");
			printf("-----------------------------------------------------------------------------------------");
			printf("\n");
			

			/*---For implicit calculations, load the previous displacement state---*/
			
			if(Precice_IsReadCheckpointRequired())
			{
    			//Precice_ReadIterationCheckpoint(&simulationData, v );
        		Precice_FulfilledReadCheckpoint();
    		}
			

    		memcpy(&sti[0],&stx[0],sizeof(double)*6*mi[0]*ne0);
		
    		++*kode;

			


		}	// Linear static loop ends here

   		/* for cyclic symmetric sectors: duplicating the results */

    	if(*mcs>0)
		{
		ptime=*ttime+time;
      	frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,t1act,
		   fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
                   nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,
                   qfn,ialset,istartset,iendset,trab,inotr,ntrans,orab,
	           ielorien,norien,sti,veold,&noddiam,set,nset,emn,thicke,
	           jobnamec,&ne0,cdn,mortar,nmat,qfx,ielprop,prop);
    	}

    	else
		{
			if(strcmp1(&filab[1044],"ZZS")==0)
			{
	    		NNEW(neigh,ITG,40**ne);
	    		NNEW(ipneigh,ITG,*nk);
			}
			ptime=*ttime+time;
			frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    	kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    	nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    	ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    	mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    	cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    	thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,ielprop,prop);

			if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
    	}	

    	/* updating the .sta file */

    	iitsta=1;
    	FORTRAN(writesta,(istep,&iinc,&icutb,&iitsta,ttime,&time,&dtime));







    	
		SFREE(v);
		SFREE(stn);
		SFREE(inum);
    	SFREE(b);
		SFREE(stx);
		SFREE(fn);

    	if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
    	if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
    	if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
    	if(strcmp1(&filab[2175],"CONT")==0) SFREE(cdn);

		

  	} // Implicit loop ends here


  	

  	/* updating the loading at the end of the step; 
     important in case the amplitude at the end of the step
     is not equal to one */

  	for(k=0;k<*nboun;++k){xbounold[k]=xbounact[k];}
  	for(k=0;k<*nforc;++k){xforcold[k]=xforcact[k];}
  	for(k=0;k<2**nload;++k){xloadold[k]=xloadact[k];}
  	for(k=0;k<7**nbody;k=k+7){xbodyold[k]=xbodyact[k];}

  	SFREE(xbounact);
	SFREE(xforcact);
	SFREE(xloadact);
	SFREE(t1act);
	SFREE(ampli);
  	SFREE(xbodyact);
	if(*nbody>0) SFREE(ipobody);SFREE(xstiff);

  	if(iglob!=0)
	{
		SFREE(integerglob);
		SFREE(doubleglob);
	}

  	*irowp=irow;*enerp=ener;*xstatep=xstate;*ipkonp=ipkon;*lakonp=lakon;
  	*konp=kon;*ielmatp=ielmat;*ielorienp=ielorien;*icolp=icol;

  	(*ttime)+=(*tper);

    /* Adapter: Free the memory */
  	Precice_FreeData( &simulationData );

	/* gkdas2: write final updated coordinates and free memory */
	printf("Writing displaced coordinates to file\n");
	writeUpdatedCO(coUpdated, *nk, mt); 
	printf("Writing baseline coordinates\n");
	writeBaselineCO(co, *nk,mt);
	SFREE(coUpdated);	
  return;
}
