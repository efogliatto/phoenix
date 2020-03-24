#include "periodicXY3D.H"


void periodicXY3D( basicMesh& mesh, const uint nx, const uint ny, const uint nz ) {

    
    // Only 2 boundaries: Z0 and Z1
    
    mesh.bd.nbd = 2;
    
    mesh.bd.nbdelem.resize(2);
    
    mesh.bd.bdPoints.resize(2);


    
    

    // Number of elements per boundary
    
    mesh.bd.nbdelem[0] = nx*ny;

    mesh.bd.nbdelem[1] = nx*ny;
        
    mesh.bd.bdPoints[0].resize( mesh.bd.nbdelem[0] );
    
    mesh.bd.bdPoints[1].resize( mesh.bd.nbdelem[1] );



    

    // Boundary elements asignment

    uint k = 0;

    for( uint j = 0 ; j < ny ; j++ ) {
    
	for( uint i = 0 ; i < nx ; i++ ) {

	    mesh.bd.bdPoints[0][k] = i + j*nx;

	    mesh.bd.bdPoints[1][k] = i + j*nx + (nz-1)*nx*ny;

	    k++;

	}

    }



    

    // Neighbour correction    


    // Move over X-boundary points

    // i = 0, i = nx-1
    // 0 < j < ny
    // 1 < k < nz - 1
    
    for( uint k = 1 ; k < (nz-1) ; k++ ) {
    
    	for( uint j = 0 ; j < ny ; j++ ) {

	    
	    uint idx_0 = j*nx + k*nx*ny;

	    uint idx_1 = idx_0 + nx - 1;


    	    for( uint vid = 0 ; vid < mesh.Q ; vid++ ) {
		
		
    		// X = 0
		
    		if(mesh.nb[idx_0][vid] == -1) {

    		    mesh.nb[idx_0][vid] = mesh.nb[idx_1][vid];

    		}


    		// X = nx

    		if(mesh.nb[idx_1][vid] == -1) {

    		    mesh.nb[idx_1][vid] = mesh.nb[idx_0][vid];

    		}
	    
    	    }

    	}

    }




    // Move over Y-boundary points

    // 0 < i < nx
    // j = 0, j = ny -1
    // 1 < k < nz - 1
    
    for( uint k = 1 ; k < (nz-1) ; k++ ) {
    
    	for( uint i = 0 ; i < nx ; i++ ) {

	    
	    uint idx_0 = i  +  k*nx*ny;

	    uint idx_1 = idx_0 + (ny-1)*nx;

	    
    	    for( uint vid = 0 ; vid < mesh.Q ; vid++ ) {

		
		
    		// Y = 0
		
    		if(mesh.nb[idx_0][vid] == -1) {

    		    mesh.nb[idx_0][vid] = mesh.nb[idx_1][vid];

    		}


    		// Y = ny

    		if(mesh.nb[idx_1][vid] == -1) {

    		    mesh.nb[idx_1][vid] = mesh.nb[idx_0][vid];

    		}
	    
    	    }

    	}

    }    





    // Corner correction on Z = 0

    uint cid[4];

    cid[0] = 0;

    cid[1] = nx-1;

    cid[2] = nx*(ny-1);

    cid[3] = nx*ny -1;


    for( uint i = 0 ; i < 4 ; i++ ) {

    	for( uint j = 0 ; j < 4 ; j++ ) {

    	    if( i!=j )  {

    		for( uint vid = 1 ; vid < mesh.Q ; vid++ ) {

    		    if (  ( mesh.nb[ cid[i] ][vid] == -1)  &&  (mesh.nb[ cid[j] ][vid] != -1)  ){

    			mesh.nb[ cid[i] ][vid] = mesh.nb[cid[j] ][vid];

    		    }
		
    		}
    	    }
	    
    	}

    }

    



    // Corner correction on Z = nz - 1 

    cid[0] = 0 + (nz-1)*nx*ny;

    cid[1] = nx-1  + (nz-1)*nx*ny;

    cid[2] = nx*(ny-1)  + (nz-1)*nx*ny;

    cid[3] = nx*ny -1  + (nz-1)*nx*ny;


    for( uint i = 0 ; i < 4 ; i++ ) {

    	for( uint j = 0 ; j < 4 ; j++ ) {

    	    if( i!=j )  {

    		for( uint vid = 1 ; vid < mesh.Q ; vid++ ) {

    		    if (  ( mesh.nb[ cid[i] ][vid] == -1)  &&  (mesh.nb[ cid[j] ][vid] != -1)  ){

    			mesh.nb[ cid[i] ][vid] = mesh.nb[cid[j] ][vid];

    		    }
		
    		}
    	    }
	    
    	}

    }


    
}
