#include <stdlib.h>
#include <stdio.h>
#include <basicMesh.h>

#ifdef USE_METIS


void kmetisDecomp( uint* owner, basicMesh* mesh, uint np )  {


    if( np > 1 ) {

	// Create graph for metis decomposition.

    
	// Count total number of edges (counted twice)

	uint nedges = 0;

	uint i, k;

	for( i = 0 ; i < mesh->nPoints ; i++ ) {

	    for( k = 1 ; k < mesh->Q ; k++ ) {

		if( mesh->nb[i][k] != -1 ) {

		    nedges++;

		}

	    }

	}


	nedges = nedges / 2;



	// Write graph file

	FILE* gfile = fopen("lattice/lattice.graph","w");

    
	fprintf(gfile,"%d %d\n", mesh->nPoints, nedges);

	for( i = 0 ; i < mesh->nPoints ; i++ ) {

	    for( k = 1 ; k < mesh->Q ; k++ ) {

		if( mesh->nb[i][k] != -1 ) {

		    fprintf(gfile,"%d ", mesh->nb[i][k]+1);

		}

	    }

	    fprintf(gfile,"\n");

	}
    
    
	fclose(gfile);

    




	// Apply metis algorithm. Call external function gpmetis

	char cmd[100];

	sprintf(cmd,"gpmetis lattice/lattice.graph %d > log.gpmetis",np);

	uint status = system( cmd );


	if (!status) {


	    // Read gpmetis result and load into owner

	    sprintf(cmd,"lattice/lattice.graph.part.%d",np);

	    gfile = fopen(cmd,"r");

	
	    for( i = 0 ; i < mesh->nPoints ; i++ ) {

		status = fscanf(gfile,"%d",&owner[i]);

	    }
	
	
	    fclose(gfile);
	

	    /* status = system("rm lattice/lattice.graph*"); */

	
	}

	else {

	    printf("\n   [ERROR]  gpmetis not executed\n\n");

	    exit(1);

	}

    
    }



    

    else {

	uint i;

	for ( i = 0 ; i < mesh->nPoints ; i++ ) {

	    owner[i] = 0;

	}

    }

}


#endif





#ifndef USE_METIS

void kmetisDecomp( uint* owner, basicMesh* mesh, uint np )  {

    printf("\n   [ERROR]  METIS library not included\n\n");

    exit(1);

}

#endif
