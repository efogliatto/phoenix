#include <latticeMesh.H>

#include <fstream>

#include <string>

#include <algorithm>


using namespace std;


/** Default constructor */

latticeMesh::latticeMesh( const int& pid, const bool& msg ) : parallel(pid) {

    if(msg) {

	if(pid==0) {

	    cout << "Reading lattice mesh" << endl << endl;

	}	

    }
    
    // Mesh points

    readPoints();


    // Lattice Neighbours

    readNeighbours();


    // Boundary nodes

    readBoundaryNodes();

    

}


/** Default destructor */

latticeMesh::~latticeMesh() {}



/** Read mesh points */

const void latticeMesh::readPoints() {


    const int pid = parallel.id();

    
    ifstream inFile( ("processor" + to_string(pid) + "/lattice/points").c_str() );

    if( inFile.is_open() ) {


	// Resize arrays
    
	nPoints = parallel.local() + parallel.ghosts();

	points.resize(nPoints);

	for(uint i = 0 ; i < nPoints ; i++)
	    points[i].resize(3);

    

	// Read all points

	int aux;

	inFile >> aux;

	for( uint i = 0 ; i < parallel.local() ; i++ ) {

	    inFile >> points[i][0];

	    inFile >> points[i][1];

	    inFile >> points[i][2];

	}


	inFile.close();	

    }

    else {

	if( pid == 0 ) {

	    cout << endl << " [ERROR] Unable to open file" << "processor" + to_string(pid) + "/lattice/points"  << endl << endl;

	}

    }
    



    
    // Read ghosts

    inFile.open( ("processor" + to_string(pid) + "/lattice/ghosts").c_str() );

    if( inFile.is_open() ) {


	int aux;

	inFile >> aux;

	for( uint i = parallel.local() ; i < nPoints ; i++ ) {

	    inFile >> points[i][0];

	    inFile >> points[i][1];

	    inFile >> points[i][2];

	}

	inFile.close();
	

    }

    else {

	if( pid == 0 ) {

	    cout << endl << " [ERROR] Unable to open file" << "processor" + to_string(pid) + "/lattice/ghosts"  << endl << endl;

	}

    }
    
    


}





/** Read lattice neighbours */

const void latticeMesh::readNeighbours() {


    const int pid = parallel.id();

    
    ifstream inFile( ("processor" + to_string(pid) + "/lattice/neighbours").c_str() );

    if( inFile.is_open() ) {

	
	// Read lattice properties (local points, d and q)
	
	uint np, d, q;

	inFile >> np;

	inFile >> d;

	inFile >> q;


	
	// Lattice model creation

	latticeModelCreator lbm;
	
	lbmodel = lbm.create(d,q);

	    

	// Resize neighbour array
	
	nb.resize(np);

	for( uint i = 0 ; i < np ; i++ )
	    nb[i].resize(q);
	

	for( uint i = 0 ; i < np ; i++ ) {

	    for( uint j = 0 ; j < q ; j++ ) {
	    
		inFile >> nb[i][j];

	    }

	}


	inFile.close();	

    }

    else {

	if( pid == 0 ) {

	    cout << endl << " [ERROR] Unable to open file" << "processor" + to_string(pid) + "/lattice/neigbours"  << endl << endl;

	}

    }
    

}






/** Read boundary nodes */

const void latticeMesh::readBoundaryNodes() {


    const int pid = parallel.id();

    
    ifstream inFile( ("processor" + to_string(pid) + "/lattice/boundary").c_str() );

    if( inFile.is_open() ) {

	// Total number of boundaries

	uint nob;

	inFile >> nob;


	// Read boundaries

	for( uint i = 0 ; i < nob ; i++ ) {


	    // Boundary name and size

	    string bdname;

	    inFile >> bdname;

	    uint bdsize;

	    inFile >> bdsize;


	    // Resize map and read

	    boundary[bdname].resize(bdsize);

	    for( uint j = 0 ; j < bdsize ; j++ )
		inFile >> boundary[bdname][j];
	    
	    
	}
	

	inFile.close();	

    }

    else {

	if( pid == 0 ) {

	    cout << endl << " [ERROR] Unable to open file" << "processor" + to_string(pid) + "/lattice/boundary"  << endl << endl;

	}

    }



    // For nodes on boundary, create map with boundary name

    for( const auto &bd : boundary ) {

	for( auto id : bd.second ) {

	    nodeToBnd[id] = bd.first;

	}

    }    



    // Check if node is on boundary

    isOnBnd.resize(local(),0);

    for(uint i = 0 ; i < local() ; i++) {

	if( nodeToBnd.find(i) != nodeToBnd.end() )
	    isOnBnd[i] = 1;

    }


    // Create map for closests nodes

    const uint q = lmodel()->q();    
    
    for( const auto &bd : boundary ) {

    	for( auto id : bd.second ) {



	    // Dont check for periodic nodes

	    bool is_periodic(true);
	    
    	    for( uint k = 0 ; k < q ; k++ ) {

		if(nb[id][k] == -1)
		    is_periodic = false;

	    }

	    

	    // Check neighbour of boundary node

	    if( !is_periodic ) {
	    
		for( uint k = 0 ; k < q ; k++ ) {

		    int nbid = nb[id][k];

		
		    if(  (nbid != -1)  &&  (nbid < (int)local())  ) {

			bool is_on_bnd(false);
		    
			for( uint l = 0 ; l < q ; l++ ) {

			    if( nb[nbid][l] == -1 ) {

				is_on_bnd = true;

			    }

			}


			// Neighbour on boundary node is not on boundary

			if( !is_on_bnd ) {


			    // Before adding to closestNodes, check if it is not already in the map

			    bool isClosest(false);
			
			    for( const auto &cnbd : closestNodes ) {

				if( std::find(cnbd.second.begin(), cnbd.second.end(), nbid) != cnbd.second.end() )
				    isClosest = true;
			    
			    }

			    if( !isClosest )
				closestNodes[bd.first].push_back(nbid);
			

			}

		    


		    }
		
		}
		
    	    }

	    

    	}

    }    
    
    

}






/** Reference to boundary indices */

const vector<uint>& latticeMesh::boundaryNodes( const std::string& bdname ) const {

    
    // Check if boundary nodes are in map

    if ( boundary.find(bdname) != boundary.end() ) {

	return boundary.at(bdname);
	
    }

    else {

	cout << " [ERROR]  Unable to find nodes from " << bdname << endl;

	exit(1);

    }

    return boundary.at("");


}




/** Boundary which node belongs */

const string latticeMesh::nodeToBoundary( const uint id ) const {

    string bd("");
	
    if( nodeToBnd.find(id) != nodeToBnd.end() ) {
	
	bd = nodeToBnd.at(id);

    }

    return bd;

}
