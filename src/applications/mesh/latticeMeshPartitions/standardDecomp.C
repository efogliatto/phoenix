#include <stdlib.h>
#include <stdio.h>
#include <basicMesh.h>

void standardDecomp( uint* owner, basicMesh* mesh, uint np )  {

    
    uint* procBegin = (uint*)malloc( np * sizeof(uint) );
    uint* procEnd   = (uint*)malloc( np * sizeof(uint) );

    uint i,
	 pid;

    
    // Split on processors
    for( i = 0 ; i < np ; i++ ) {
	
    	if(i != 0) {
	    
    	    procBegin[i] = procEnd[i-1] + 1;
    	    procEnd[i] = procBegin[i] + (int)(mesh->nPoints/np) - 1;
	    
    	}
	
    	else {
	    
    	    procBegin[i] = 0;
    	    procEnd[0] = (int)(mesh->nPoints/np) - 1;
	    
    	}
	
    }
    
    procEnd[np - 1] = mesh->nPoints - 1;

    

    
    // Ownership vector
    
    for( pid = 0 ; pid < np ; pid++ ) {
	
    	for( i = procBegin[pid] ; i <= procEnd[pid] ; i++ ) {
	    
    	    owner[i] = pid;
	    
    	}

    }   

    
}
