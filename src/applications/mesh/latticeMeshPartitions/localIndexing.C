// #include <stdlib.h>
// #include <stdio.h>
// #include <basicMesh.h>
// #include <basic.h>

#include "localIndexing.H"

using namespace std;


void localIndexing ( basicMesh& mesh, vector< vector<int> >& local, vector< vector<int> >& nGhosts, vector<uint>& owner, uint np ) {

    
    
    // Create an auxiliary copy of local indexing (in order to track changes)
    
    vector<int> aux( mesh.nPoints );
    

    
    // Check all points for every processor
    
    for( uint j = 0 ; j < np ; j++ ) {

	
    	uint id = 0;


    	// Assign local nodes first
	
    	for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	    if( owner[i] == j ) {

    		local[i][j] = id;

    		id++;

    	    }

    	}



    	// Total number of local elements
	
    	nGhosts[j][0] = id;
	


    	// Copy to aux
	
    	for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	    aux[i] = local[i][j];

    	}


	

    	// Check between neighbours using only local points

    	for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	    if(  aux[i] != -1 ) {
		
    		for( uint velId = 0 ; velId < mesh.Q ; velId++ ) {
		    
    		    int nid = mesh.nb[i][velId];
		    
    		    if( nid != -1 )  {
			
    			if( local[nid][j] == -1 ) {

    			    local[nid][j] =  id;

    			    id++;

    			}

    		    }

    		}

    	    }

    	}



    	// Total number of ghost elements
	
    	nGhosts[j][1] = id - nGhosts[j][0];
	

    }
    
    
    

}
