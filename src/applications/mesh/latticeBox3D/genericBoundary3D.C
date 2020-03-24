#include "genericBoundary3D.H"


void genericBoundary3D( basicMesh& mesh, const uint nx, const uint ny, const uint nz ) {
    
    mesh.bd.nbd = 6;
    
    mesh.bd.nbdelem.resize(6);
    mesh.bd.bdPoints.resize(6);

    mesh.bd.nbdelem[0] = ny*nz;
    mesh.bd.nbdelem[1] = ny*nz;
    mesh.bd.nbdelem[2] = (nx-2)*nz;
    mesh.bd.nbdelem[3] = (nx-2)*nz;
    mesh.bd.nbdelem[4] = (nx-2)*(ny-2);
    mesh.bd.nbdelem[5] = (nx-2)*(ny-2);

    
    mesh.bd.bdPoints[0].resize( mesh.bd.nbdelem[0] );
    mesh.bd.bdPoints[1].resize( mesh.bd.nbdelem[1] );
    mesh.bd.bdPoints[2].resize( mesh.bd.nbdelem[2] );
    mesh.bd.bdPoints[3].resize( mesh.bd.nbdelem[3] );   
    mesh.bd.bdPoints[4].resize( mesh.bd.nbdelem[4] );
    mesh.bd.bdPoints[5].resize( mesh.bd.nbdelem[5] );   


    



    // Move over X-boundary points

    int id = 0;
    
    for( uint k = 0 ; k < nz ; k++ ) {

	for( uint j = 0 ; j < ny ; j++ ) {

	    mesh.bd.bdPoints[0][id] = j*nx + k*nx*ny;
	
	    mesh.bd.bdPoints[1][id] = (nx-1) + j*nx + k*nx*ny;

	    id++;

	}
	
    }
    


    
    // Move over Y-boundary points

    id = 0;
    
    for( uint k = 0 ; k < nz ; k++ ) {
    
	for( uint i = 1 ; i < (nx-1) ; i++ ) {

	    mesh.bd.bdPoints[2][id] = i + k*nx*ny; 
	
	    mesh.bd.bdPoints[3][id] = i + (ny-1)*nx + k*nx*ny;

	    id++;

	}

    }

    

    // Move over Z-boundary points

    id = 0;
    
    for( uint j = 1 ; j < (ny-1) ; j++ ) {
    
    	for( uint i = 1 ; i < (nx-1) ; i++ ) {

    	    mesh.bd.bdPoints[4][id] = i + j*nx;
	
    	    mesh.bd.bdPoints[5][id] = i + j*nx + (nz-1)*nx*ny;

	    id++;
	    
    	}

    }


    
}
