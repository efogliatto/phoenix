#include <stdio.h>
#include <stdlib.h>
#include <basicMesh.h>


void periodicY( basicMesh* mesh, uint nx, uint ny ) {


    mesh->bd.nbd = 4;
    
    mesh->bd.nbdelem = (uint*)malloc( 4 * sizeof(uint) );
    mesh->bd.bdPoints = (uint**)malloc( 4 * sizeof(uint*) );

    mesh->bd.nbdelem[0] = (ny-2);
    mesh->bd.nbdelem[1] = (ny-2);
    mesh->bd.nbdelem[2] = nx;
    mesh->bd.nbdelem[3] = nx;
    
    mesh->bd.bdPoints[0] = (uint*)malloc( mesh->bd.nbdelem[0] * sizeof(uint) );
    mesh->bd.bdPoints[1] = (uint*)malloc( mesh->bd.nbdelem[1] * sizeof(uint) );
    mesh->bd.bdPoints[2] = (uint*)malloc( mesh->bd.nbdelem[2] * sizeof(uint) );
    mesh->bd.bdPoints[3] = (uint*)malloc( mesh->bd.nbdelem[3] * sizeof(uint) );   


    
    int j,velId;


    // Move over X-boundary points
    for( j = 1 ; j < (ny-1) ; j++ ) {

    	mesh->bd.bdPoints[0][j-1] = j*nx;
	
    	mesh->bd.bdPoints[1][j-1] = (nx-1) + j*nx;

    }
    

    
    // Move over Y-boundary points
    for( j = 0 ; j < nx ; j++ ) {

    	mesh->bd.bdPoints[2][j] = j;
	
    	mesh->bd.bdPoints[3][j] = j + (ny-1)*nx;


    	// Assign periodic neighbours
    	for( velId = 0 ; velId < mesh->Q ; velId++ ) {

    	    if(mesh->nb[j][velId] == -1) {

    		mesh->nb[j][velId] = mesh->nb[ j + (nx-1)*ny ][velId];

    	    }

    	    if(mesh->nb[ j + (nx-1)*ny ][velId] == -1) {

    		mesh->nb[ j + (nx-1)*ny ][velId] = mesh->nb[j][velId];

    	    }
	    
    	}

    }

    

    
    
}
