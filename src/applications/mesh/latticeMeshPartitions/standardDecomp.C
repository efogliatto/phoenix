#include <iostream>

#include "standardDecomp.H"

#include <vector>

#include <fstream>

using namespace std;


void standardDecomp( vector<uint>& owner, basicMesh& mesh, uint np )  {

    
    vector<uint> procBegin( np );

    vector<uint> procEnd( np );

    
    // Split on processors
    for( uint i = 0 ; i < np ; i++ ) {
	
    	if(i != 0) {
	    
    	    procBegin[i] = procEnd[i-1] + 1;
    	    procEnd[i] = procBegin[i] + (int)(mesh.nPoints/np) - 1;
	    
    	}
	
    	else {
	    
    	    procBegin[i] = 0;
    	    procEnd[0] = (int)(mesh.nPoints/np) - 1;
	    
    	}
	
    }
    
    procEnd[np - 1] = mesh.nPoints - 1;

    

    
    // Ownership vector
    
    for( uint pid = 0 ; pid < np ; pid++ ) {
	
    	for( uint i = procBegin[pid] ; i <= procEnd[pid] ; i++ ) {
	    
    	    owner[i] = pid;
	    
    	}

    }

    

    // Write lattice.graph.part for other purposes

    ofstream outFile;

    outFile.open( "lattice/lattice.graph.part." + std::to_string(np) );


    for( const auto o : owner )
	outFile << o << endl;

    
}
