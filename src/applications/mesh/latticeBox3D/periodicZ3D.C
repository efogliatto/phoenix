#include "periodicZ3D.H"


void periodicZ3D( basicMesh& mesh, const uint nx, const uint ny, const uint nz ) {

    
    // Only 4 boundaries: X0, X1, Y0 and Y1

    mesh.bd.nbd = 4;
    
    mesh.bd.nbdelem.resize(4);
    
    mesh.bd.bdPoints.resize(4);


    
    

    // Number of elements per boundary
    
    mesh.bd.nbdelem[0] = nx*ny;

    mesh.bd.nbdelem[1] = nx*ny;

    mesh.bd.nbdelem[2] = (nx-2)*ny;

    mesh.bd.nbdelem[3] = (nx-2)*ny;

    for( uint i = 0 ; i < mesh.bd.nbd ; i++ )    
	mesh.bd.bdPoints[i].resize( mesh.bd.nbdelem[i] );

   




    // Boundary elements asignment

    // X-boundary

    uint id = 0;

    for( uint k = 0 ; k < nz ; k++ ) {
    
	for( uint j = 0 ; j < ny ; j++ ) {

	    mesh.bd.bdPoints[0][id] = j*nx + k*nx*ny;

	    mesh.bd.bdPoints[1][id] = (nx-1) + j*nx + k*nx*ny;

	    id++;

	}

    }



    // Y-boundary

    id = 0;

    for( uint k = 0 ; k < nz ; k++ ) {
    
	for( uint i = 1 ; i < (nx-1) ; i++ ) {

	    mesh.bd.bdPoints[2][id] = i + k*nx*ny;

	    mesh.bd.bdPoints[3][id] = i  +  (ny-1)*nx  +  k*nx*ny;

	    id++;

	}

    }    

    



    // Neighbour correction    


    // Move over Z-boundary points

    // 1 < i < nx - 1
    // 1 < j < ny - 1
    // k = 0, k = nz - 1
    
    
    for( uint j = 1 ; j < (ny-1) ; j++ ) {
    
    	for( uint i = 1 ; i < (nx-1) ; i++ ) {

	    
	    uint idx_0 = i + j*nx;

	    uint idx_1 = idx_0 + (nz-1)*nx*ny;



    	    for( uint vid = 0 ; vid < mesh.Q ; vid++ ) {
		
		
    		// Z = 0
		
    		if(mesh.nb[idx_0][vid] == -1) {

    		    mesh.nb[idx_0][vid] = mesh.nb[idx_1][vid];

    		}


    		// Z = nz - 1

    		if(mesh.nb[idx_1][vid] == -1) {

    		    mesh.nb[idx_1][vid] = mesh.nb[idx_0][vid];

    		}
	    
    	    }

    	}

    }
    




    // Correct periodicity on y edges in z-direction

    for( uint i = 0 ; i < nx ; i++ ) {

	uint idx_0 = i;

	uint idx_1 = idx_0 + (nz-1)*nx*ny;

	for( uint vid = 0 ; vid < mesh.Q ; vid++ ) {
		
	    if(mesh.nb[idx_0][vid] == -1) {

		mesh.nb[idx_0][vid] = mesh.nb[idx_1][vid];

	    }

	    if(mesh.nb[idx_1][vid] == -1) {

		mesh.nb[idx_1][vid] = mesh.nb[idx_0][vid];

	    }
	    
	}



	idx_0 = i + (ny-1)*nx;

	idx_1 = idx_0 + (nz-1)*nx*ny;

	for( uint vid = 0 ; vid < mesh.Q ; vid++ ) {
		
	    if(mesh.nb[idx_0][vid] == -1) {

		mesh.nb[idx_0][vid] = mesh.nb[idx_1][vid];

	    }

	    if(mesh.nb[idx_1][vid] == -1) {

		mesh.nb[idx_1][vid] = mesh.nb[idx_0][vid];

	    }
	    
	}
	

    }










    // Correct periodicity on x edges in z-direction

    for( uint j = 1 ; j < (ny-1) ; j++ ) {

	uint idx_0 = j*ny;

	uint idx_1 = idx_0 + (nz-1)*nx*ny;

	for( uint vid = 0 ; vid < mesh.Q ; vid++ ) {
		
	    if(mesh.nb[idx_0][vid] == -1) {

		mesh.nb[idx_0][vid] = mesh.nb[idx_1][vid];

	    }

	    if(mesh.nb[idx_1][vid] == -1) {

		mesh.nb[idx_1][vid] = mesh.nb[idx_0][vid];

	    }
	    
	}



	idx_0 = (nx-1) + j*nx;

	idx_1 = idx_0 + (nz-1)*nx*ny;

	for( uint vid = 0 ; vid < mesh.Q ; vid++ ) {
		
	    if(mesh.nb[idx_0][vid] == -1) {

		mesh.nb[idx_0][vid] = mesh.nb[idx_1][vid];

	    }

	    if(mesh.nb[idx_1][vid] == -1) {

		mesh.nb[idx_1][vid] = mesh.nb[idx_0][vid];

	    }
	    
	}
	

    }

    
}
