#include "genericBoundary.H"


void genericBoundary( basicMesh& mesh, const uint nx, const uint ny ) {
    
    mesh.bd.nbd = 4;
    
    mesh.bd.nbdelem.resize(4);
    mesh.bd.bdPoints.resize(4);

    mesh.bd.nbdelem[0] = ny;
    mesh.bd.nbdelem[1] = ny;
    mesh.bd.nbdelem[2] = (nx-2);
    mesh.bd.nbdelem[3] = (nx-2);

    mesh.bd.bdPoints[0].resize( mesh.bd.nbdelem[0] );
    mesh.bd.bdPoints[1].resize( mesh.bd.nbdelem[1] );
    mesh.bd.bdPoints[2].resize( mesh.bd.nbdelem[2] );
    mesh.bd.bdPoints[3].resize( mesh.bd.nbdelem[3] );   


    

    // Move over X-boundary points
    for( uint j = 0 ; j < ny ; j++ ) {

    	mesh.bd.bdPoints[0][j] = j*nx;
	
    	mesh.bd.bdPoints[1][j] = (nx-1) + j*nx;

    }
    
    
    // Move over Y-boundary points
    for( uint j = 1 ; j < (nx-1) ; j++ ) {

    	mesh.bd.bdPoints[2][j-1] = j;
	
    	mesh.bd.bdPoints[3][j-1] = j + (nx-1)*ny;

    }


    
}
