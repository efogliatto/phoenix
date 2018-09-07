#include <mpiInfo.H>

#include <fstream>

#include <sstream>

#include <string>

#include <iostream>


using namespace std;


/** Default constructor */

mpiInfo::mpiInfo( const uint& id ) : pid(id) {


    // Read recv ghosts
    
    ifstream inFile( ("processor" + to_string(pid) + "/lattice/recvGhosts").c_str() );

    if( inFile.is_open() ) {

	
	// Total number of processors

	inFile >> worldSize;


	// Resize identifiers

	recvGhosts.resize(worldSize);


	// Read elements

	for( uint i = 0 ; i < worldSize ; i++ ) {

	    uint gid, nrg;

	    inFile >> gid;
	    
	    inFile >> nrg;

	    recvGhosts[gid].resize(nrg);

	    for(uint j = 0 ; j < nrg ; j++)
		inFile >> recvGhosts[gid][j];

	}
 

	inFile.close();

    }

    else {

	if( pid == 0 ) {

	    cout << endl << " [ERROR] Unable to open file" << "processor" + to_string(pid) + "/lattice/recvGhosts"  << endl << endl;

	}

    }







    // Read send ghosts
    
    inFile.open( ("processor" + to_string(pid) + "/lattice/sendGhosts").c_str() );

    if( inFile.is_open() ) {

	
	// Total number of processors

	inFile >> worldSize;


	// Resize identifiers

	sendGhosts.resize(worldSize);


	// Read elements

	for( uint i = 0 ; i < worldSize ; i++ ) {

	    uint gid, nsg;

	    inFile >> gid;

	    inFile >> nsg;

	    sendGhosts[gid].resize(nsg);

	    for(uint j = 0 ; j < nsg ; j++)
		inFile >> sendGhosts[gid][j];

	}
 

	inFile.close();

    }

    else {

	if( pid == 0 ) {

	    cout << endl << " [ERROR] Unable to open file" << "processor" + to_string(pid) + "/lattice/sendGhosts"  << endl << endl;

	}

    }    



    

    // Total number of elements per patch

    nodesPerPatch.resize( worldSize );

    for( uint i = 0 ; i < worldSize ; i++ ) {	

	inFile.open( ("processor" + to_string(pid) + "/lattice/points").c_str() );

	uint npp;

	inFile >> npp;

	inFile.close();

	nodesPerPatch[i] = npp;

	if( i == pid ) {

	    nlocal = npp;
	    
	}

	

	inFile.open( ("processor" + to_string(pid) + "/lattice/ghosts").c_str() );

	inFile >> npp;

	inFile.close();

	nodesPerPatch[i] += npp;


	if( i == pid ) {

	    nghosts = npp;
	    
	}

    }







    
    

}



/** Default destructor */

mpiInfo::~mpiInfo() {}
