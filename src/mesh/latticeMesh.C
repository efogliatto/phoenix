#include <latticeMesh.H>

#include <fstream>

#include <string>


using namespace std;


/** Default constructor */

latticeMesh::latticeMesh( const int& pid ) : parallel(pid) {

    // Mesh points

    readPoints();

    

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
