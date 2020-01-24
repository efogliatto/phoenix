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

	uint clevel( (uint)pdict.lookUpOrDefault<scalar>("coarseLevel",np) );


       
	
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


	idx_t options[METIS_NOPTIONS];

	METIS_SetDefaultOptions(options);       

	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

	options[METIS_OPTION_NUMBERING] = 0;

	options[METIS_OPTION_CONTIG] = 1;

	options[METIS_OPTION_UFACTOR] = 1;




	

	int status =  METIS_PartMeshNodal( &ne, &nn, eptr, eind, NULL, NULL, &nparts, NULL, options, &objval, epart, npart);



	if(status == METIS_OK) {
	    
	    cout << "Finished METIS nodal decomposition" << endl << endl;

	}

	else {

    	    cout << "\n   [ERROR]  Unsuccesful METIS decomposition\n\n";

    	    exit(1);

	}






	// Partition clustering:
	// partitions are reorderer in clevel groups

	// First set communication matrix

	vector< vector<uint> > commMat( np );

	for( uint i = 0 ; i < np ; i++ )
	    commMat[i].resize(np,0);


	// Move over points and add weight

	for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

	    for( uint k = 1 ; k < mesh.Q ; k++ ) {

		int nbid = mesh.nb[i][k];

		if( nbid != -1 )
		    commMat[ npart[i] ][ npart[nbid] ]++;

	    }

	}




	// Create graph for partition communication

	vector<uint> newOwnerId(np);
	

	if( clevel > 1 ) {

	    // Count adjacency

	    uint npadj(0);

	    for( uint i = 0 ; i < np ; i++ ) {

		for( uint j = 0 ; j < np ; j++ ) {

		    if( i != j ) {

			if( commMat[i][j] != 0 )
			    npadj++;

		    }

		}

	    }



	    // Allocate graph info

	    idx_t* xadj = (idx_t*)malloc( (np+1) * sizeof(idx_t) );

	    idx_t* adjncy = (idx_t*)malloc( npadj * sizeof(idx_t) );

	    uint count(0);

	    for( uint i = 0 ; i < np ; i++ ) {

		xadj[i] = count;

		for( uint j = 0 ; j < np ; j++ ) {

		    if( i != j ) {

			if( commMat[i][j] != 0 ) {

			    adjncy[count] = j;

			    count++;
			    
			}

		    }

		}

	    }

	    xadj[np] = npadj;




	    
	    
	    // Apply graph decomposition

	    idx_t nvtxs = np;

	    idx_t ncon = 1;

	    idx_t nparts_proc = clevel;

	    idx_t objval_proc;

	    idx_t* part_proc = (idx_t*)malloc( np * sizeof(idx_t) );		

	    status =  METIS_PartGraphRecursive( &nvtxs, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts_proc, NULL, NULL, NULL, &objval_proc, part_proc );




	    if(status == METIS_OK) {
	    
	    	cout << "Finished graph partitioning for decomposed domain" << endl << endl;

	    }

	    else {

	    	cout << "\n   [ERROR]  Unsuccesful METIS decomposition\n\n";

	    	exit(1);

	    }


	    


	    // Reassign processor number

	    {

		uint procId(0);

		for( uint i = 0 ; i < clevel ; i++ ) {

		    for( uint j = 0 ; j < np ; j++ ) {

			if( (uint)part_proc[j] == i ) {

			    newOwnerId[j] = procId;

			    procId++;

			}

		    }

		}

		    
	    }



	    // Release memory

	    free(xadj);
	    free(adjncy);
	    free(part_proc);

    

	}


       
	

	


	// Copy partition to owner array

	if( clevel > 1 ) {

	    for( uint i = 0 ; i < mesh.nPoints ; i++ )	    
		owner[i] = newOwnerId[ npart[i] ];

	}

	else {

	    for( uint i = 0 ; i < mesh.nPoints ; i++ )	    
		owner[i] = npart[i];

	}
	

	// Deallocate ugly memory

	free(eptr);
	free(eind);
	free(epart);
	free(npart);




    
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
