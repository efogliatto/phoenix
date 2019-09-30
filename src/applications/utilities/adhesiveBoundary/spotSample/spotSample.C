#include <spotSample.H>

#include <dictionary.H>


using namespace std;


/** Default constructor */

spotSample::spotSample(const std::string& bdname) {


    // Read node indices for boundary "bdname"

    ifstream inFile;

    inFile.open( "lattice/boundary");

    if( inFile.is_open() == false ){
	
    	cout << " [ERROR]  Unable to find file lattice/boundary" << endl;
	
    	exit(1);
	
    }


    // Total number of boundaries
    
    uint nbds;

    inFile >> nbds;

    
    for(uint bnd = 0 ; bnd < nbds ; bnd++) {

	string auxName;

	uint nnodes;

	uint auxnodes;
	

	inFile >> auxName;

	inFile >> nnodes;


	if( auxName == bdname ) {

	    _nodes.resize(nnodes);
	    
	    for(uint j = 0 ; j < nnodes ; j++) {
	    
		inFile >> _nodes[j];

	    }

	}

	else {

	    for(uint j = 0 ; j < nnodes ; j++) {
	    
		inFile >> auxnodes;

	    }

	}
	

    }


    inFile.close();





    // Look for nodes location

    _location.resize( _nodes.size() );

    for(uint i = 0 ; i < _location.size() ; i++)
	_location[i].resize(3);
    

    inFile.open( "lattice/points");

    if( inFile.is_open() == false ){
	
    	cout << " [ERROR]  Unable to find file lattice/points" << endl;
	
    	exit(1);
	
    }


    uint npoints;

    inFile >> npoints;

    vector< vector<uint> > points;

    points.resize(npoints);

    for(uint i = 0 ; i < npoints ; i++)
	points[i].resize(3);


    for(uint i = 0 ; i < npoints ; i++) {
	
    	inFile >> points[i][0];

    	inFile >> points[i][1];

    	inFile >> points[i][2];	

    }



    // Match node location

    for(uint i = 0 ; i < _nodes.size() ; i++) {

	for(uint j = 0 ; j < 3 ; j++)
	    _location[i][j] = points[ _nodes[i] ][j];

    }
    

}


/** Destructor */

spotSample::~spotSample() {}
