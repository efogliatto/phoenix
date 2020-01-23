#include <fstream>

#include <iostream>

#include "mpmetisDecomp.H"


// #ifdef USE_METIS

using namespace std;

void mpmetisDecomp( vector<uint>& owner, basicMesh& mesh, uint np )  {


    if( np > 1 ) {
	

    	// // Create graph for metis decomposition.

    
    	// // Count total number of edges (counted twice)

    	// uint nedges = 0;

    	// for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	//     for( uint k = 1 ; k < maxneigh ; k++ ) {

    	// 	if( mesh.nb[i][k] != -1 ) {

    	// 	    nedges++;

    	// 	}

    	//     }

    	// }


    	// nedges = nedges / 2;



    // 	// Write graph file

    // 	ofstream gfile;

    // 	gfile.open( "lattice/lattice.graph" );	

    
    // 	gfile << mesh.nPoints << " " << nedges << endl;

    // 	for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    // 	    for( uint k = 1 ; k < maxneigh ; k++ ) {

    // 		if( mesh.nb[i][k] != -1 ) {

    // 		    gfile << mesh.nb[i][k]+1 << " ";

    // 		}

    // 	    }

    // 	    gfile << endl;

    // 	}
    
    
    // 	gfile.close();

    




    // 	// Apply metis algorithm. Call external function gpmetis

    // 	char cmd[100];

    // 	sprintf(cmd,"gpmetis lattice/lattice.graph %d -contig > log.gpmetis",np);

    // 	uint status = system( cmd );


    // 	if (!status) {


    // 	    // Read gpmetis result and load into owner

    // 	    sprintf(cmd,"lattice/lattice.graph.part.%d",np);

    // 	    ifstream inFile;

    // 	    inFile.open( cmd );	
	    

	
    // 	    for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    // 		inFile >> owner[i];

    // 	    }
	
	
    // 	    inFile.close();	

	
    // 	}

    // 	else {

    // 	    cout << "\n   [ERROR]  gpmetis not executed\n\n";

    // 	    exit(1);

    // 	}

    
    }



    else {

    	for ( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	    owner[i] = 0;

    	}

    }
	

}


// #endif





// #ifndef USE_METIS

// using namespace std;

// #include <iostream>

// void kmetisDecomp( vector<uint>& owner, basicMesh& mesh, uint np )  {

//     cout << "\n   [ERROR]  METIS library not included\n" << endl;

//     exit(1);

// }

// #endif
