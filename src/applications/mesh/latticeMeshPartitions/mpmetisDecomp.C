#include <fstream>

#include <iostream>

#include "mpmetisDecomp.H"

#include <dictionary.H>



#ifdef USE_METIS

#include <metis.h>


using namespace std;

void mpmetisDecomp( vector<uint>& owner, basicMesh& mesh, uint np )  {


    if( np > 1 ) {
	

	// Read multi-level information

	dictionary pdict( "properties/parallel" );

	uint coarseLevel( (uint)pdict.lookUpOrDefault<scalar>("coarseLevel",np) );


       
	
    	// Create mesh for metis decomposition.

    
    	// Resize METIS mesh arrays

    	uint nadj = mesh.ncells * mesh.cellType;

	idx_t* eptr = (idx_t*)malloc( (mesh.ncells+1) * sizeof(idx_t) );

	idx_t* eind = (idx_t*)malloc( nadj * sizeof(idx_t) );


	// Fill METIS mesh arrays

	for( uint i = 0 ; i < mesh.ncells ; i++ ) {

	    eptr[i] = (idx_t)i*mesh.cellType;


	    for( uint j = 0 ; j < mesh.cellType ; j++ )
	    	eind[ eptr[i] + j ] = mesh.vtkCells[i][j];

	}

	eptr[mesh.ncells] = nadj;



	

	// Apply METIS algorithm

	idx_t ne = mesh.ncells;

	idx_t nn = mesh.nPoints;

	idx_t nparts = np;

	idx_t objval;

	idx_t* epart = (idx_t*)malloc( ne * sizeof(idx_t) );

	idx_t* npart = (idx_t*)malloc( nn * sizeof(idx_t) );		

	int status =  METIS_PartMeshNodal( &ne, &nn, eptr, eind, NULL, NULL, &nparts, NULL, NULL, &objval, epart, npart);



	if(status) {
	    
	    cout << "Finished METIS nodal decomposition" << endl << endl;

	}

	else {

    	    cout << "\n   [ERROR]  Unsuccesful METIS decomposition\n\n";

    	    exit(1);

	}






	// Copy partition to owner array

	for( uint i = 0 ; i < mesh.nPoints ; i++ )	    
	    owner[i] = npart[i];



    




    // 	// Apply metis algorithm. Call external function gpmetis

    // 	char cmd[100];

    // 	sprintf(cmd,"gpmetis lattice/lattice.graph %d -contig > log.gpmetis",np);

    // 	uint status = system( cmd );


    // 	if (!status) {


    // 	    // Read gpmetis result and load into owner

    // 	    sprintf(cmd,"lattice/lattice.graph.part.%d",np);

    // 	    ifstream inFile;

    // 	    inFile.open( cmd );	
	    

	
    // 	    for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    // 		inFile >> owner[i];

    // 	    }
	
	
    // 	    inFile.close();	

	
    // 	}

    // 	else {

    // 	    cout << "\n   [ERROR]  gpmetis not executed\n\n";

    // 	    exit(1);

    // 	}

    
    }



    else {

    	for ( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	    owner[i] = 0;

    	}

    }
	

}


#endif





#ifndef USE_METIS

using namespace std;

#include <iostream>

void kmetisDecomp( vector<uint>& owner, basicMesh& mesh, uint np )  {

    cout << "\n   [ERROR]  METIS library not included\n" << endl;

    exit(1);

}


#endif
