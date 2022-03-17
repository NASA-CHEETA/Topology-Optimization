/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#include <stdlib.h>
#include <assert.h>
#include "PreciceInterface.h"
#include "ConfigReader.h"
#include "precice/SolverInterfaceC.h"

void Precice_Setup( char * configFilename, char * participantName, SimulationData * sim )
{
  assert(sim != NULL);
  assert(configFilename != NULL);
  assert(participantName != NULL);

	//printf( "Setting up preCICE participant %s, using config file: %s\n", participantName, configFilename );
	fflush( stdout );

	// Read the YAML config file
	
  AdapterConfig adapterConfig;
	ConfigReader_Read( configFilename, participantName, &adapterConfig);
	

  assert(adapterConfig.interfaces != NULL);
  assert(adapterConfig.preciceConfigFilename != NULL);
  assert(adapterConfig.numInterfaces > 0);

  sim->numPreciceInterfaces = adapterConfig.numInterfaces;

	// Create the solver interface and configure it - Alex: Calculix is always a serial participant (MPI size 1, rank 0)
	precicec_createSolverInterface( participantName, adapterConfig.preciceConfigFilename, 0, 1 );

	// Create interfaces as specified in the config file
	sim->preciceInterfaces = (struct PreciceInterface**) calloc( adapterConfig.numInterfaces, sizeof( PreciceInterface* ) );

	for(int i = 0 ; i < adapterConfig.numInterfaces; i++ )
	{
    	InterfaceConfig * config = adapterConfig.interfaces + i;
		
		sim->preciceInterfaces[i] = malloc( sizeof( PreciceInterface ) );
		
		PreciceInterface_Create( sim->preciceInterfaces[i], sim, config );
		
	}

  // At this point we are done with the configuration
  AdapterConfig_Free(&adapterConfig);

  

	// Initialize variables needed for the coupling
	NNEW( sim->coupling_init_v, double, sim->mt * sim->nk );

	// Initialize preCICE
	sim->precice_dt = precicec_initialize();

	// Initialize coupling data
	Precice_InitializeData( sim );
}

void Precice_InitializeData( SimulationData * sim )
{
	printf( "Initializing coupling data\n" );
	fflush( stdout );

	Precice_WriteCouplingData( sim );
	precicec_initialize_data();
	Precice_ReadCouplingData( sim );
}

void Precice_AdjustSolverTimestep( SimulationData * sim )
{
	if( isSteadyStateSimulation( sim->nmethod ) )
	{
		printf("Adjusting time step for linear static analysis \n" );
		fflush( stdout );

		// For steady-state simulations, we will always compute the converged steady-state solution in one coupling step
		/*
		*sim->theta = 0;
		*sim->tper = 1;
		*sim->dtheta = 1;
		*/
		// Set the solver time step to be the same as the coupling time step
		printf("Setting solver_dt (syncronization interval) to %f \n", sim->precice_dt);
		sim->solver_dt = sim->precice_dt;
	}
	else
	{
		/*---Uncomment below for a pure dynamic solution---*/
		/*
		printf( "Adjusting time step for transient step\n" );
		printf( "precice_dt dtheta = %f, dtheta = %f, solver_dt = %f\n", sim->precice_dt / *sim->tper, *sim->dtheta, fmin( sim->precice_dt, *sim->dtheta * *sim->tper ) );
		fflush( stdout );

		// Compute the normalized time step used by CalculiX
		*sim->dtheta = fmin( sim->precice_dt / *sim->tper, *sim->dtheta );

		// Compute the non-normalized time step used by preCICE
		sim->solver_dt = ( *sim->dtheta ) * ( *sim->tper );

		*/

		/*---Uncomment below for a pure steady state solution---*/
		printf( "Adjusting time step for linear static analysis\n" );
		fflush(stdout);

		/*---For steady-state simulations, we will always compute the converged steady-state solution in one coupling step---*/
		*sim->theta = 0;      /*---Set sum of all previous increments to zero---*/
		*sim->tper = 1;		  /*---Set step size to 1---*/
		*sim->dtheta = 1;	  /*---Set increment size to one---*/

		/*---Set the solver time step to be the same as the coupling time step---*/
	    	
		//sim->solver_dt = 0;
		sim->solver_dt = sim->precice_dt;





	}
}

void Precice_Advance( SimulationData * sim )
{
	printf( "Adapter calling advance()...\n" );
	fflush( stdout );
    
	printf("Advancing with solver dt: %f \n", sim->solver_dt);
	sim->precice_dt = precicec_advance(1);
	//sim->precice_dt = precicec_advance( sim->solver_dt );
}

bool Precice_IsCouplingOngoing()
{
	return precicec_isCouplingOngoing();
}

bool Precice_IsReadCheckpointRequired()
{
	return precicec_isActionRequired( "read-iteration-checkpoint" );
}

bool Precice_IsWriteCheckpointRequired()
{
	return precicec_isActionRequired( "write-iteration-checkpoint" );
}

void Precice_FulfilledReadCheckpoint()
{
	precicec_markActionFulfilled( "read-iteration-checkpoint" );
}

void Precice_FulfilledWriteCheckpoint()
{
	precicec_markActionFulfilled( "write-iteration-checkpoint" );
}

void Precice_ReadIterationCheckpoint( SimulationData * sim, double * v )
{

	printf( "Adapter reading checkpoint...\n" );
	fflush( stdout );

	// Reload time
  //	*( sim->theta ) = sim->coupling_init_theta;

	// Reload step size
//	*( sim->dtheta ) = sim->coupling_init_dtheta;

	// Reload solution vector v
	memcpy( v, sim->coupling_init_v, sizeof( double ) * sim->mt * sim->nk );
}

void Precice_WriteIterationCheckpoint( SimulationData * sim, double * v )
{

	printf( "Adapter writing checkpoint...\n" );
	fflush( stdout );

	// Save time
	//sim->coupling_init_theta = *( sim->theta );

	// Save step size
//	sim->coupling_init_dtheta = *( sim->dtheta );

	// Save solution vector v
	memcpy( sim->coupling_init_v, v, sizeof( double ) * sim->mt * sim->nk );
}

void Precice_ReadCouplingData( SimulationData * sim )
{

	printf( "Adapter reading coupling data...\n" );
	
	fflush( stdout );


    
	PreciceInterface ** interfaces = sim->preciceInterfaces;

	printf("Failing here \n");
	
	int numInterfaces = sim->numPreciceInterfaces;
	int i, j;
	
  // Extract force Data



	if( precicec_isReadDataAvailable() )
	{
		for( i = 0 ; i < numInterfaces ; i++ )
		{

			for( j = 0 ; j < interfaces[i]->numReadData ; j++ )
			{

				switch( interfaces[i]->readData[j] )
				{
				case TEMPERATURE:
					// Read and set temperature BC
					precicec_readBlockScalarData( interfaces[i]->temperatureDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeScalarData );
					setNodeTemperatures( interfaces[i]->nodeScalarData, interfaces[i]->numNodes, interfaces[i]->xbounIndices, sim->xboun );
					printf( "Reading TEMPERATURE coupling data with ID '%d'. \n",interfaces[i]->temperatureDataID );
					break;
				case HEAT_FLUX:
					// Read and set heat flux BC
					precicec_readBlockScalarData( interfaces[i]->fluxDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData );
					setFaceFluxes( interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload );
					printf( "Reading HEAT_FLUX coupling data with ID '%d'. \n",interfaces[i]->fluxDataID );
					break;
				case SINK_TEMPERATURE:
					// Read and set sink temperature in convective film BC
					precicec_readBlockScalarData( interfaces[i]->kDeltaTemperatureReadDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData );
					setFaceSinkTemperatures( interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload );
					printf( "Reading SINK_TEMPERATURE coupling data with ID '%d'. \n",interfaces[i]->kDeltaTemperatureReadDataID );
					break;
				case HEAT_TRANSFER_COEFF:
					// Read and set heat transfer coefficient in convective film BC
					precicec_readBlockScalarData( interfaces[i]->kDeltaReadDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData );
					setFaceHeatTransferCoefficients( interfaces[i]->faceCenterData, interfaces[i]->numElements, interfaces[i]->xloadIndices, sim->xload );
					printf( "Reading HEAT_TRANSFER_COEFF coupling data with ID '%d'. \n",interfaces[i]->kDeltaReadDataID );
					break;
        		case FORCES:
					// Read and set forces as concentrated loads (Neumann BC)
					precicec_readBlockVectorData( interfaces[i]->forcesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
					setNodeForces( interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->dim, interfaces[i]->xforcIndices, sim->xforc);
					printf( "Reading FORCES coupling data with ID '%d'. \n",interfaces[i]->forcesDataID );
          			int Number_nodes = interfaces[i]->numNodes;
          			printf( "Number of Interface nodes '%d'. \n", Number_nodes);
        //    for ( int k = 0; k < 3 * Number_nodes; k++)
        //   {
      //      printf( "Nodal frc '%lf'. \n", interfaces[i]->nodeVectorData[k]);
    //       }
          double *Fx;
          double *Fy;
          double *Fz;
          double *NID;
          double *F;
          F = (double*) malloc(3*Number_nodes*sizeof(double));
          NID = (int*) malloc(Number_nodes*sizeof(double));
          Fx = (double*) malloc(Number_nodes*sizeof(double));
          Fy = (double*) malloc(Number_nodes*sizeof(double));
          Fz = (double*) malloc(Number_nodes*sizeof(double));

          for ( int k = 0; k < 3 * Number_nodes; k++)
          {
            F[k] =  interfaces[i]->nodeVectorData[k];
          }

          for ( int k = 0; k < Number_nodes; k++)
          {
            NID[k] =  interfaces[i]->nodeIDs[k];
          }

          // Assign the first and last values for Fx

          Fx[0] = F[0];
          Fy[0] = F[1];
          Fz[0] = F[2];

          Fx[Number_nodes-1] = F[3*Number_nodes-3];
          Fy[Number_nodes-1] = F[3*Number_nodes-2];
          Fz[Number_nodes-1] = F[3*Number_nodes-1];

		  double sum_x = 0.0;
		  double sum_y = 0.0;
		  double sum_z = 0.0;

		  int Count_x  = 0;
		  int Count_y  = 0;
		  int Count_z  = 0;
		  double Fx_Avg, Fy_Avg, Fz_Avg;

          for (int n = 1; n <Number_nodes-2; n++)
          {
              Fx[n] = F[3*n];
			  sum_x +=Fx[n];
			  Count_x++;
          }

		  Fx_Avg = (sum_x/Count_x);


          for (int n = 1; n < Number_nodes-2; n++)
          {
            Fy[n] = F[3*n+1];
			sum_y +=Fy[n];
			Count_y++;
          }
		  Fy_Avg = (sum_y/Count_y);

          for (int n = 1; n < Number_nodes-2; n++)
          {
            Fz[n] = F[3*n+2];
			sum_z +=Fz[n];
			Count_z++;
          }
		  Fz_Avg = (sum_z/Count_z);

		  /*---Print average Forces in all dimensions---*/

          printf("\n");
		  printf("---------------------------------------------------------------------------\n");
		  printf("Average aerodynamic traction in x: %f \n", Fx_Avg);
		  printf("Average aerodynamic traction in y: %f \n", Fy_Avg);
		  printf("Average aerodynamic traction in z: %f \n", Fz_Avg);
		  printf("---------------------------------------------------------------------------\n");
		  printf("\n");

          FILE *fptr;
          fptr = fopen("Forces.csv", "w");
          if (fptr == NULL)
          {
            printf("ERROR! Unable to open write-file.\n");
            exit(1);
          }
            fprintf(fptr, " Surface Node ID  Fx (N)  Fy (N) Fz (N) \n");
          for ( int aa = 0; aa < Number_nodes ; aa ++)
          {
            fprintf(fptr, " %lf , %lf , %lf  , %lf \n", NID[aa], Fx[aa], Fy[aa], Fz[aa]);
          }
           fclose(fptr);


        //  for ( int k = 0; k < Number_nodes; k++)
        //  {
        //    printf(  " %lf, %lf, %lf, %lf \n", NID[k], Fx[k], Fy[k], Fz[k]);
        //  }
          free(NID);
          free(Fx);
          free(Fy);
          free(Fz);
          free(F);

          //  for ( int k = 0; k < Number_nodes; k++)
          // {
          //  printf( "Nodal ID Global'%d'. \n", interfaces[i]->nodeIDs[k]);
        //   }
					break;
				case DISPLACEMENTS:
					// Read and set displacements as single point constraints (Dirichlet BC)
					precicec_readBlockVectorData( interfaces[i]->displacementsDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
					setNodeDisplacements( interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->dim, interfaces[i]->xbounIndices, sim->xboun );
					printf( "Reading DISPLACEMENTS coupling data with ID '%d'. \n",interfaces[i]->displacementsDataID );
					break;
				case DISPLACEMENTDELTAS:
					printf( "DisplacementDeltas cannot be used as read data\n" );
					fflush( stdout );
					exit( EXIT_FAILURE );
					break;
				case VELOCITIES:
					printf( "Velocities cannot be used as read data\n" );
					fflush( stdout );
					exit( EXIT_FAILURE );
					break;
				case POSITIONS:
					printf( "Positions cannot be used as read data.\n" );
					fflush( stdout );
					exit( EXIT_FAILURE );
					break;
				}
			}
		}
	}
}





void Precice_WriteCouplingData( SimulationData * sim )
{

	printf( "Adapter writing coupling data...\n" );
	fflush( stdout );

	PreciceInterface ** interfaces = sim->preciceInterfaces;
	int numInterfaces = sim->numPreciceInterfaces;
	int i, j;
	int iset;

	if( precicec_isWriteDataRequired( sim->solver_dt ) || precicec_isActionRequired( "write-initial-data" ) )
	{
		for( i = 0 ; i < numInterfaces ; i++ )
		{

			for( j = 0 ; j < interfaces[i]->numWriteData ; j++ )
			{

				switch( interfaces[i]->writeData[j] )
				{
				case TEMPERATURE:
					getNodeTemperatures( interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->mt, interfaces[i]->nodeScalarData );
					precicec_writeBlockScalarData( interfaces[i]->temperatureDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeScalarData );
					printf( "Writing TEMPERATURE coupling data with ID '%d'. \n",interfaces[i]->temperatureDataID );
					break;
				case HEAT_FLUX:
					iset = interfaces[i]->faceSetID + 1; // Adjust index before calling Fortran function
					FORTRAN( getflux, ( sim->co,
										sim->ntmat_,
										sim->vold,
										sim->cocon,
										sim->ncocon,
										&iset,
										sim->istartset,
										sim->iendset,
										sim->ipkon,
										*sim->lakon,
										sim->kon,
										sim->ialset,
										sim->ielmat,
										sim->mi,
										interfaces[i]->faceCenterData
										)
							 );
					precicec_writeBlockScalarData( interfaces[i]->fluxDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, interfaces[i]->faceCenterData );
					printf( "Writing HEAT_FLUX coupling data with ID '%d'. \n",interfaces[i]->fluxDataID );
					break;
				case SINK_TEMPERATURE:
					iset = interfaces[i]->faceSetID + 1; // Adjust index before calling Fortran function
					double * myKDelta = malloc( interfaces[i]->numElements * sizeof( double ) );
					double * T = malloc( interfaces[i]->numElements * sizeof( double ) );
					FORTRAN( getkdeltatemp, ( sim->co,
											  sim->ntmat_,
											  sim->vold,
											  sim->cocon,
											  sim->ncocon,
											  &iset,
											  sim->istartset,
											  sim->iendset,
											  sim->ipkon,
											  *sim->lakon,
											  sim->kon,
											  sim->ialset,
											  sim->ielmat,
											  sim->mi,
											  myKDelta,
											  T
											  )
							 );
					precicec_writeBlockScalarData( interfaces[i]->kDeltaTemperatureWriteDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, T );
					printf( "Writing SINK_TEMPERATURE coupling data with ID '%d'. \n",interfaces[i]->kDeltaTemperatureWriteDataID );
					free( T );
					break;
				case HEAT_TRANSFER_COEFF:
					precicec_writeBlockScalarData( interfaces[i]->kDeltaWriteDataID, interfaces[i]->numElements, interfaces[i]->preciceFaceCenterIDs, myKDelta );
					printf( "Writing HEAT_TRANSFER_COEFF coupling data with ID '%d'. \n",interfaces[i]->kDeltaWriteDataID );
					free( myKDelta );
					break;
				case DISPLACEMENTS:
					getNodeDisplacements( interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->vold, sim->mt, interfaces[i]->nodeVectorData );
					precicec_writeBlockVectorData( interfaces[i]->displacementsDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
					printf( "Writing DISPLACEMENTS coupling data with ID '%d'. \n",interfaces[i]->displacementsDataID );
					break;
				case DISPLACEMENTDELTAS:
					getNodeDisplacementDeltas( interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->vold, sim->coupling_init_v, sim->mt, interfaces[i]->nodeVectorData );
					precicec_writeBlockVectorData( interfaces[i]->displacementDeltasDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
					printf( "Writing DISPLACEMENTDELTAS coupling data with ID '%d'. \n",interfaces[i]->displacementDeltasDataID );
					break;
				case VELOCITIES:
					getNodeVelocities( interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->veold, sim->mt, interfaces[i]->nodeVectorData );
					precicec_writeBlockVectorData( interfaces[i]->velocitiesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
					printf( "Writing VELOCITIES coupling data with ID '%d'. \n",interfaces[i]->velocitiesDataID );
					break;
				case POSITIONS:
					getNodeCoordinates( interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->co, sim->vold, sim->mt, interfaces[i]->nodeVectorData );
					precicec_writeBlockVectorData( interfaces[i]->positionsDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
					printf( "Writing POSITIONS coupling data with ID '%d'. \n",interfaces[i]->positionsDataID );
					break;
				case FORCES:
					getNodeForces( interfaces[i]->nodeIDs, interfaces[i]->numNodes, interfaces[i]->dim, sim->fn, sim->mt, interfaces[i]->nodeVectorData );
					precicec_writeBlockVectorData( interfaces[i]->forcesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
					printf( "Writing FORCES coupling data with ID '%d'. \n",interfaces[i]->forcesDataID );
					break;
				}
			}
		}

		if( precicec_isActionRequired( "write-initial-data" ) )
		{
			precicec_markActionFulfilled( "write-initial-data" );
		}
	}
}




void Precice_FreeData( SimulationData * sim )
{
	int i;

	if( sim->coupling_init_v != NULL ){
		free( sim->coupling_init_v );
	}

	for( i = 0 ; i < sim->numPreciceInterfaces ; i++ )
	{
		PreciceInterface_FreeData( sim->preciceInterfaces[i] );
		if( sim->preciceInterfaces[i] != NULL ){
			free( sim->preciceInterfaces[i] );
		}
	}

	precicec_finalize();
}

void PreciceInterface_Create( PreciceInterface * interface, SimulationData * sim, InterfaceConfig const * config )
{
	
	interface->dim = precicec_getDimensions();

	

	// Initialize pointers as NULL
	interface->elementIDs = NULL;
	interface->faceIDs = NULL;
	interface->faceCenterCoordinates = NULL;
	interface->preciceFaceCenterIDs = NULL;
	interface->nodeCoordinates = NULL;
	interface->preciceNodeIDs = NULL;
	interface->triangles = NULL;
	interface->nodeScalarData = NULL;
	interface->nodeVectorData = NULL;
	interface->faceCenterData = NULL;
	interface->xbounIndices = NULL;
	interface->xloadIndices = NULL;
	interface->xforcIndices = NULL;

	

  // Initialize data ids to -1
	interface->temperatureDataID = -1;
	interface->fluxDataID = -1;
	interface->kDeltaWriteDataID = -1;
	interface->kDeltaTemperatureWriteDataID = -1;
	interface->kDeltaReadDataID = -1;
	interface->kDeltaTemperatureReadDataID = -1;
	interface->displacementsDataID = -1;
	interface->displacementDeltasDataID = -1;
	interface->positionsDataID = -1;
	interface->velocitiesDataID = -1;
	interface->forcesDataID = -1;

	

	//Mapping Type
	// The patch identifies the set used as interface in Calculix
	interface->name = strdup( config->patchName );
	// Calculix needs to know if nearest-projection mapping is implemented. config->map = 1 is for nearest-projection, config->map = 0 is for everything else
	interface->mapNPType = config->map;
    
	// Nodes mesh
	interface->nodesMeshID = -1;
	interface->nodesMeshName = NULL;
  if ( config->nodesMeshName ) 
  {
    interface->nodesMeshName = strdup( config->nodesMeshName );
    PreciceInterface_ConfigureNodesMesh( interface, sim );
  }

  

	// Face centers mesh
	interface->faceCentersMeshID = -1;
	interface->faceCentersMeshName = NULL;

  if ( config->facesMeshName ) 
  {
	  interface->faceCentersMeshName = strdup( config->facesMeshName );
		//Only configure a face center mesh if necesary; i.e. do not configure it for FSI simulations, also do not configure tetra faces if no face center mesh is used (as in FSI simulations)
    PreciceInterface_ConfigureFaceCentersMesh( interface, sim );
		// Triangles of the nodes mesh (needs to be called after the face centers mesh is configured!)
    PreciceInterface_ConfigureTetraFaces( interface, sim );
  }

  

	PreciceInterface_ConfigureCouplingData( interface, sim, config );
	
}

void PreciceInterface_ConfigureFaceCentersMesh( PreciceInterface * interface, SimulationData * sim )
{
	//printf("Entering ConfigureFaceCentersMesh \n");
	char * faceSetName = toFaceSetName( interface->name );
	interface->faceSetID = getSetID( faceSetName, sim->set, sim->nset );
	interface->numElements = getNumSetElements( interface->faceSetID, sim->istartset, sim->iendset );

	interface->elementIDs = malloc( interface->numElements * sizeof( ITG ) );
	interface->faceIDs = malloc( interface->numElements * sizeof( ITG ) );
	getSurfaceElementsAndFaces( interface->faceSetID, sim->ialset, sim->istartset, sim->iendset, interface->elementIDs, interface->faceIDs );

	interface->faceCenterCoordinates = malloc( interface->numElements * 3 * sizeof( double ) );
	interface->preciceFaceCenterIDs = malloc( interface->numElements * 3 * sizeof( int ) );
	getTetraFaceCenters( interface->elementIDs, interface->faceIDs, interface->numElements, sim->kon, sim->ipkon, sim->co, interface->faceCenterCoordinates, interface->preciceFaceCenterIDs );


	interface->faceCentersMeshID = precicec_getMeshID( interface->faceCentersMeshName );
	interface->preciceFaceCenterIDs = malloc( interface->numElements * sizeof( int ) );

	precicec_setMeshVertices( interface->faceCentersMeshID, interface->numElements, interface->faceCenterCoordinates, interface->preciceFaceCenterIDs);

}

void PreciceInterface_ConfigureNodesMesh( PreciceInterface * interface, SimulationData * sim )
{

	//printf("Entering configureNodesMesh \n");
	char * nodeSetName = toNodeSetName( interface->name );
	interface->nodeSetID = getSetID( nodeSetName, sim->set, sim->nset );
	interface->numNodes = getNumSetElements( interface->nodeSetID, sim->istartset, sim->iendset );
	//printf("numNodes = %d \n", interface->numNodes);
	interface->nodeIDs = &sim->ialset[sim->istartset[interface->nodeSetID] - 1]; //Lucia: make a copy

	interface->nodeCoordinates = malloc( interface->numNodes * interface->dim * sizeof( double ) );
	getNodeCoordinates( interface->nodeIDs, interface->numNodes, interface->dim, sim->co, sim->vold, sim->mt, interface->nodeCoordinates );

	if( interface->nodesMeshName != NULL )
	{
		//printf("nodesMeshName is not null \n");
		interface->nodesMeshID = precicec_getMeshID( interface->nodesMeshName );
		interface->preciceNodeIDs = malloc( interface->numNodes * sizeof( int ) );
		//interface->preciceNodeIDs = malloc( interface->numNodes * 3 * sizeof( int ) );
		//getNodeCoordinates( interface->nodeIDs, interface->numNodes, sim->co, sim->vold, sim->mt, interface->nodeCoordinates, interface->preciceNodeIDs );
		precicec_setMeshVertices( interface->nodesMeshID, interface->numNodes, interface->nodeCoordinates, interface->preciceNodeIDs );
	}

	if (interface->mapNPType == 1)
	{
			PreciceInterface_NodeConnectivity( interface, sim );
	}
}

void PreciceInterface_NodeConnectivity( PreciceInterface * interface, SimulationData * sim )
{
	int numElements;
	char * faceSetName = toFaceSetName( interface->name );
	interface->faceSetID = getSetID( faceSetName, sim->set, sim->nset );
	numElements = getNumSetElements( interface->faceSetID, sim->istartset, sim->iendset );
	interface->triangles = malloc( numElements * 3 * sizeof( ITG ) );
	interface->elementIDs = malloc( numElements * sizeof( ITG ) );
	interface->faceIDs = malloc( numElements * sizeof( ITG ) );
	interface->faceCenterCoordinates = malloc( numElements * 3 * sizeof( double ) );
	getSurfaceElementsAndFaces( interface->faceSetID, sim->ialset, sim->istartset, sim->iendset, interface->elementIDs, interface->faceIDs );
	interface->numElements = numElements;
	interface->triangles = malloc( numElements * 3 * sizeof( ITG ) );
	PreciceInterface_ConfigureTetraFaces( interface, sim );
}

void PreciceInterface_EnsureValidNodesMeshID( PreciceInterface * interface )
{
	if( interface->nodesMeshID < 0 )
	{
		printf( "Nodes mesh not provided in YAML config file\n" );
		fflush( stdout );
		exit( EXIT_FAILURE );
	}
}

void PreciceInterface_ConfigureTetraFaces( PreciceInterface * interface, SimulationData * sim )
{
	int i;
	printf("Setting node connectivity for nearest projection mapping: \n");
	if( interface->nodesMeshName != NULL )
	{
		interface->triangles = malloc( interface->numElements * 3 * sizeof( ITG ) );
		getTetraFaceNodes( interface->elementIDs, interface->faceIDs,  interface->nodeIDs, interface->numElements, interface->numNodes, sim->kon, sim->ipkon, interface->triangles );

		for( i = 0 ; i < interface->numElements ; i++ )
		{
			precicec_setMeshTriangleWithEdges( interface->nodesMeshID, interface->triangles[3*i], interface->triangles[3*i+1], interface->triangles[3*i+2] );
		}
	}
}

void PreciceInterface_ConfigureCouplingData( PreciceInterface * interface, SimulationData * sim, InterfaceConfig const * config )
{

	interface->nodeScalarData = malloc( interface->numNodes * sizeof( double ) );
	interface->nodeVectorData = malloc( interface->numNodes * 3 * sizeof( double ) );
	interface->faceCenterData = malloc( interface->numElements * sizeof( double ) );

	int i;

        interface->numReadData = config->numReadData;
		
        if (config->numReadData > 0) interface->readData = malloc( config->numReadData * sizeof( int ) );
		
	for(i = 0 ; i < config->numReadData ; i++ )
	{

		if( isEqual( config->readDataNames[i], "Temperature" ) )
		{
			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->readData[i] = TEMPERATURE;
			interface->xbounIndices = malloc( interface->numNodes * sizeof( int ) );
			interface->temperatureDataID = precicec_getDataID( "Temperature", interface->nodesMeshID );
			getXbounIndices( interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices, TEMPERATURE );
			printf( "Read data '%s' found with ID # '%d'.\n", config->readDataNames[i], interface->temperatureDataID);
		}
		else if ( isEqual( config->readDataNames[i], "Heat-Flux" ) )
		{
			interface->readData[i] = HEAT_FLUX;
			interface->xloadIndices = malloc( interface->numElements * sizeof( int ) );
			getXloadIndices( "DFLUX", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices );
			interface->fluxDataID = precicec_getDataID( "Heat-Flux", interface->faceCentersMeshID );
			printf( "Read data '%s' found with ID # '%d'.\n", config->readDataNames[i],interface->fluxDataID );
		}
		else if ( startsWith( config->readDataNames[i], "Sink-Temperature-" ) )
		{
			interface->readData[i] = SINK_TEMPERATURE;
			interface->xloadIndices = malloc( interface->numElements * sizeof( int ) );
			getXloadIndices( "FILM", interface->elementIDs, interface->faceIDs, interface->numElements, sim->nload, sim->nelemload, sim->sideload, interface->xloadIndices );
			interface->kDeltaTemperatureReadDataID = precicec_getDataID( config->readDataNames[i], interface->faceCentersMeshID );
			printf( "Read data '%s' found with ID # '%d'.\n", config->readDataNames[i], interface->kDeltaTemperatureReadDataID);
		}
		else if ( startsWith( config->readDataNames[i], "Heat-Transfer-Coefficient-" ) )
		{
			interface->readData[i] = HEAT_TRANSFER_COEFF;
			interface->kDeltaReadDataID = precicec_getDataID( config->readDataNames[i], interface->faceCentersMeshID );
			printf( "Read data '%s' found with ID # '%d'.\n", config->readDataNames[i],interface->kDeltaReadDataID );
		}
		else if ( startsWith( config->readDataNames[i], "Forces" ) )
		{
			
			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->readData[i] = FORCES;
			
			interface->xforcIndices = malloc( interface->numNodes * 3 * sizeof( int ) );
			
			interface->forcesDataID = precicec_getDataID( config->readDataNames[i], interface->nodesMeshID );
			
			getXforcIndices( interface->nodeIDs, interface->numNodes, sim->nforc, sim->ikforc, sim->ilforc, interface->xforcIndices );

			printf( "Read data '%s' found with ID # '%d'.\n", config->readDataNames[i],interface->forcesDataID );
		}
		else if ( startsWith( config->readDataNames[i], "Displacements" ) )
		{
			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->readData[i] = DISPLACEMENTS;
			interface->xbounIndices = malloc( interface->numNodes * 3 * sizeof( int ) );
			interface->displacementsDataID = precicec_getDataID( config->readDataNames[i], interface->nodesMeshID );
			getXbounIndices( interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices, DISPLACEMENTS );
			printf( "Read data '%s' found with ID # '%d'.\n", config->readDataNames[i], interface->displacementsDataID );
		}
		else
		{
			printf( "ERROR: Read data '%s' does not exist!\n", config->readDataNames[i] );
			exit( EXIT_FAILURE );
		}
	}

        interface->numWriteData = config->numWriteData;
        if (config->numWriteData > 0) interface->writeData = malloc( config->numWriteData * sizeof( int ) );
	for( i = 0 ; i < config->numWriteData ; i++ )
	{
		if( isEqual( config->writeDataNames[i], "Temperature" ) )
		{
			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->writeData[i] = TEMPERATURE;
			interface->temperatureDataID = precicec_getDataID( "Temperature", interface->nodesMeshID );
			printf( "Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i],interface->temperatureDataID );
		}
		else if ( isEqual( config->writeDataNames[i], "Heat-Flux" ) )
		{
			interface->writeData[i] = HEAT_FLUX;
			interface->fluxDataID = precicec_getDataID( "Heat-Flux", interface->faceCentersMeshID );
			printf( "Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i],interface->fluxDataID );
		}
		else if ( startsWith( config->writeDataNames[i], "Sink-Temperature-" ) )
		{
			interface->writeData[i] = SINK_TEMPERATURE;
			interface->kDeltaTemperatureWriteDataID = precicec_getDataID( config->writeDataNames[i], interface->faceCentersMeshID );
			printf( "Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i],interface->kDeltaTemperatureWriteDataID );
		}
		else if ( startsWith( config->writeDataNames[i], "Heat-Transfer-Coefficient-" ) )
		{
			interface->writeData[i] = HEAT_TRANSFER_COEFF;
			interface->kDeltaWriteDataID = precicec_getDataID( config->writeDataNames[i], interface->faceCentersMeshID );
			printf( "Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i],interface->kDeltaWriteDataID );
		}
		else if ( startsWith( config->writeDataNames[i], "Displacements" ) )
		{
			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->writeData[i] = DISPLACEMENTS;
			interface->displacementsDataID = precicec_getDataID( config->writeDataNames[i], interface->nodesMeshID );
			printf( "Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i],interface->displacementsDataID );
		}
		else if ( startsWith( config->writeDataNames[i], "DisplacementDeltas" ) )
		{
			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->writeData[i] = DISPLACEMENTDELTAS;
			interface->displacementDeltasDataID = precicec_getDataID( config->writeDataNames[i], interface->nodesMeshID );
			printf( "Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i],interface->displacementDeltasDataID );
		}
		else if ( startsWith( config->writeDataNames[i], "Positions" ) )
		{
			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->writeData[i] = POSITIONS;
			interface->positionsDataID = precicec_getDataID( config->writeDataNames[i], interface->nodesMeshID );
			printf( "Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i],interface->positionsDataID );
		}
		else if ( startsWith( config->writeDataNames[i], "Velocities" ) )
		{
			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->writeData[i] = VELOCITIES;
			interface->velocitiesDataID = precicec_getDataID( config->writeDataNames[i], interface->nodesMeshID );
			printf( "Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i],interface->velocitiesDataID );
		}
		else if ( startsWith( config->writeDataNames[i], "Forces" ) )
		{
			PreciceInterface_EnsureValidNodesMeshID( interface );
			interface->writeData[i] = FORCES;
			interface->forcesDataID = precicec_getDataID( config->writeDataNames[i], interface->nodesMeshID );
			printf( "Write data '%s' found with ID # '%d'.\n", config->writeDataNames[i],interface->forcesDataID );
		}
		else
		{
			printf( "ERROR: Write data '%s' does not exist!\n", config->writeDataNames[i] );
			exit( EXIT_FAILURE );
		}
	}
}

void PreciceInterface_FreeData( PreciceInterface * preciceInterface )
{
	if( preciceInterface->readData != NULL ){
		free( preciceInterface->readData );
	}

	if( preciceInterface->writeData != NULL ){
		free( preciceInterface->writeData );
	}

	if( preciceInterface->elementIDs != NULL ){
		free( preciceInterface->elementIDs );
	}

	if( preciceInterface->faceIDs != NULL ){
		free( preciceInterface->faceIDs );
	}

	if( preciceInterface->faceCenterCoordinates != NULL ){
		free( preciceInterface->faceCenterCoordinates );
	}

	if( preciceInterface->preciceFaceCenterIDs != NULL ){
		free( preciceInterface->preciceFaceCenterIDs );
	}

	if( preciceInterface->nodeCoordinates != NULL ){
		free( preciceInterface->nodeCoordinates );
	}

	if( preciceInterface->preciceNodeIDs != NULL ){
		free( preciceInterface->preciceNodeIDs );
	}

	if( preciceInterface->triangles != NULL ){
		free( preciceInterface->triangles );
	}

	if( preciceInterface->nodeScalarData != NULL ){
		free( preciceInterface->nodeScalarData );
	}

	if( preciceInterface->nodeVectorData != NULL ){
		free( preciceInterface->nodeVectorData );
	}

	if( preciceInterface->faceCenterData != NULL ){
		free( preciceInterface->faceCenterData );
	}

	if( preciceInterface->xbounIndices != NULL ){
		free( preciceInterface->xbounIndices );
	}

	if( preciceInterface->xloadIndices != NULL ){
		free( preciceInterface->xloadIndices );
	}

	if ( preciceInterface->xforcIndices != NULL ){
		free( preciceInterface->xforcIndices );
	}
}
