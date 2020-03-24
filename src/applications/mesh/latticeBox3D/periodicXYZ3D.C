#include "periodicXYZ3D.H"


void periodicXYZ3D( basicMesh& mesh, const uint nx, const uint ny, const uint nz ) {

    // No boundaries
    
    mesh.bd.nbd = 0;


    

    // Move over X-boundary points

    // i = 0, i = nx-1
    // 0 < j < ny
    // 1 < k < nz - 1
    
    for( uint k = 0 ; k < nz ; k++ ) {
    
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
    
    for( uint k = 0 ; k < nz ; k++ ) {
    
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




    // Move over Z-boundary points

    // 0 < i < nx
    // j = 0, j = ny -1
    // 1 < k < nz - 1
    
    for( uint j = 0 ; j < ny ; j++ ) {
    
    	for( uint i = 0 ; i < nx ; i++ ) {

	    
	    uint idx_0 = i  +  j*nx;

	    uint idx_1 = idx_0 + (nz-1)*nx*ny;

	    
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
    
    
}
