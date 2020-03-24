#include "periodicXY.H"


void periodicXY( basicMesh& mesh, const uint nx, const uint ny ) {


    mesh.bd.nbd = 0;
    



    // Move over X-boundary points
    for( uint j = 1 ; j < (ny-1) ; j++ ) {

    	// Assign periodic neighbours
    	for( uint velId = 0 ; velId < mesh.Q ; velId++ ) {

    	    if(mesh.nb[j*nx][velId] == -1) {

    		mesh.nb[j*nx][velId] = mesh.nb[ nx-1+j*nx  ][velId];

    	    }

    	    if(mesh.nb[nx-1+j*nx][velId] == -1) {

    		mesh.nb[nx-1+j*nx][velId] = mesh.nb[j*nx][velId];

    	    }
	    
    	}	

    }
    



    // Move over Y-boundary points
    for( uint j = 1 ; j < (nx-1) ; j++ ) {

    	// Assign periodic neighbours
    	for( uint velId = 0 ; velId < mesh.Q ; velId++ ) {

    	    if(mesh.nb[j][velId] == -1) {

    		mesh.nb[j][velId] = mesh.nb[ j + (ny-1)*nx ][velId];

    	    }

    	    if(mesh.nb[ j + (ny-1)*nx ][velId] == -1) {

    		mesh.nb[ j + (ny-1)*nx ][velId] = mesh.nb[j][velId];

    	    }
	    
    	}	

    }    



    

    // Resolve for corners
   
    int cid[4];
    cid[0] = 0;
    cid[1] = nx-1;
    cid[2] = nx*(ny-1);
    cid[3] = nx*ny -1;

    for( uint i = 0 ; i < 4 ; i++ ) {

    	for( uint j = 0 ; j < 4 ; j++ ) {

    	    if( i!=j )  {

    		for( uint velId = 1 ; velId < mesh.Q ; velId++ ) {

    		    if (  ( mesh.nb[ cid[i] ][velId] == -1)  &&  (mesh.nb[ cid[j] ][velId] != -1)  ){

    			mesh.nb[ cid[i] ][velId] = mesh.nb[cid[j] ][velId];

    		    }
		
    		}
    	    }
	    
    	}

    }

    

    
}
