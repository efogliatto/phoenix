#include "periodicX.H"


void periodicX( basicMesh& mesh, const uint nx, const uint ny ) {


    mesh.bd.nbd = 4;
    
    mesh.bd.nbdelem.resize(4);
    mesh.bd.bdPoints.resize(4);

    mesh.bd.nbdelem[0] = (ny-2);
    mesh.bd.nbdelem[1] = (ny-2);
    mesh.bd.nbdelem[2] = nx;
    mesh.bd.nbdelem[3] = nx;

    mesh.bd.bdPoints[0].resize( mesh.bd.nbdelem[0] );
    mesh.bd.bdPoints[1].resize( mesh.bd.nbdelem[1] );
    mesh.bd.bdPoints[2].resize( mesh.bd.nbdelem[2] );
    mesh.bd.bdPoints[3].resize( mesh.bd.nbdelem[3] );   


    

    // Move over X-boundary points

    for( uint j = 1 ; j < (ny-1) ; j++ ) {
	

    	mesh.bd.bdPoints[0][j-1] = j*nx;
	
    	mesh.bd.bdPoints[1][j-1] = (nx-1) + j*nx;

	
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
    
    for( uint j = 0 ; j < nx ; j++ ) {

    	mesh.bd.bdPoints[2][j] = j;
	
    	mesh.bd.bdPoints[3][j] = j + (ny-1)*nx;

    }



    
    // Corner correction

    // X0,Y0 --> X1,Y0
    
    for( uint velId = 0 ; velId < mesh.Q ; velId++ ) {

	if(mesh.nb[0][velId] == -1) {

	    mesh.nb[0][velId] = mesh.nb[(nx-1)][velId];

	}

	if(mesh.nb[(nx-1)][velId] == -1) {

	    mesh.nb[(nx-1)][velId] = mesh.nb[0][velId];

	}	

    }

    

    // X0,Y1 --> X1,Y1
    
    for( uint velId = 0 ; velId < mesh.Q ; velId++ ) {

	if(mesh.nb[nx*(ny-1)][velId] == -1) {

	    mesh.nb[nx*(ny-1)][velId] = mesh.nb[nx*ny-1][velId];

	}

	if(mesh.nb[nx*ny-1][velId] == -1) {

	    mesh.nb[nx*ny-1][velId] = mesh.nb[nx*(ny-1)][velId];

	}	

    }


    
    
}
