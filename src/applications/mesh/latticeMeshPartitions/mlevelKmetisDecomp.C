#include <fstream>

#include <iostream>

#include "mlevelKmetisDecomp.H"

#include "kmetisDecomp.H"

#include "standardDecomp.H"

#include <dictionary.H>

#include <map>


using namespace std;

void mlevelKmetisDecomp( vector<uint>& owner, basicMesh& mesh, uint np )  {


    if( np > 1 ) {


	// Read multi-level information

	dictionary pdict( "properties/parallel" );

	uint coarseLevel( (uint)pdict.lookUpOrDefault<scalar>("coarseLevel",np) );

	uint fineLevel( (uint)pdict.lookUpOrDefault<scalar>("fineLevel",1) );

	string coarseMethod( pdict.lookUpOrDefault<string>("coarseMethod","standard") );


	cout << "Multi-level decomposition" << endl;
	cout << coarseLevel << " coarse grains with " << coarseMethod << endl;
	cout << fineLevel << " fine grains with kmetis" << endl << endl;



	// Restrict neighbouring to main directions

	uint Q(5);	

	if( mesh.D == 3 )	    
	    Q = 7;





	// Decompose domain

	if( ( coarseLevel * fineLevel ) == np  ) {
	    

	    // Coarse level decomposition

	    cout << "Coarse-level decomposition" << endl << endl;
	    
	    
	    vector<uint> coarseOwner( owner.size() );

	    if( coarseMethod == "standard" ) {
		
		standardDecomp( coarseOwner, mesh, coarseLevel );

	    }

	    else {
		
		kmetisDecomp( coarseOwner, mesh, coarseLevel, Q );

	    }

	    

	    

	    // Move over coarse level and decompose into fine grains

	    cout << "Fine-level decomposition" << endl;	    

	    for( uint coarse = 0 ; coarse < coarseLevel ; coarse++ ) {


	    	cout << "  coarse level " << coarse << endl; 
		

		// Reverse indices: owner to nodes

		vector<uint> ownerToNodes;

		map<uint,uint> localMap;

		{
		    
		    uint count(0);

		    for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

			if( coarseOwner[i] == coarse ) {
			
			    ownerToNodes.push_back(i);

			    localMap[i] = count;

			    count++;

			}

		    }

		}

		


	    	// Create graph for metis decomposition.

    
	    	// Count total number of edges (counted twice)

	    	uint nedges = 0;

	    	for( const auto id : ownerToNodes ) {
		    
		    
		    for( uint k = 1 ; k < Q ; k++ ) {

			if( mesh.nb[id][k] != -1 ) {

			    if( coarseOwner[ mesh.nb[id][k] ] == coarse ) {

				nedges++;

			    }

			}

		    }

	    	}


	    	nedges = nedges / 2;





		

	    	// Write graph file

	    	ofstream gfile;

	    	gfile.open( "lattice/lattice.graph." + std::to_string(coarse) );	

    
	    	gfile << ownerToNodes.size() << " " << nedges << endl;

	    	for( const auto id : ownerToNodes ) {
		    
		    
		    for( uint k = 1 ; k < Q ; k++ ) {

			if( mesh.nb[id][k] != -1 ) {

			    if( coarseOwner[ mesh.nb[id][k] ] == coarse ) {

				gfile << localMap[ mesh.nb[id][k] ] + 1 << " ";

			    }

			}

		    }

		    gfile << endl;		

	    	}
    
    
	    	gfile.close();





	    	// Apply metis algorithm. Call external function gpmetis

	    	char cmd[100];

	    	sprintf(cmd,"gpmetis lattice/lattice.graph.%d %d > log.gpmetis",coarse,fineLevel);

	    	uint status = system( cmd );


	    	if (!status) {


	    	    // Read gpmetis result and load into owner

	    	    sprintf(cmd,"lattice/lattice.graph.%d.part.%d",coarse,fineLevel);

	    	    ifstream inFile;

	    	    inFile.open( cmd );


	
	    	    for( const auto id : ownerToNodes )			
	    		inFile >> owner[id];


	
	
	    	    inFile.close();	

	
	    	}

	    	else {

	    	    cout << "\n   [ERROR]  gpmetis not executed\n\n";

	    	    exit(1);

	    	}		


		
		

	    }

	    







	    // Finally assign global processor number

	    for( uint i = 0 ; i < mesh.nPoints ; i++ )		
		owner[i] = owner[i] + coarseOwner[i]*fineLevel;





	    

		

	    

	    

	}


	else {

	    cout << "\n\n  [ERROR]  Level mismatch " << "\n\n\n";

	    exit(1);	    

	}

	
	


    
    }



    else {

    	for ( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	    owner[i] = 0;

    	}

    }
	

}
