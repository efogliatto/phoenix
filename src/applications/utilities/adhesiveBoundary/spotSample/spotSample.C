#include <spotSample.H>

#include <dictionary.H>

#include <spotRadiusCreator.H>


using namespace std;




/** Default constructor */

spotSample::spotSample(const std::string& bdname) {


    // Nodes ids and locations

    readNodes(bdname);
    
    readLocations(bdname);


    // Assign spots positions and radius

    createSpots(bdname);


    // Create cavity model

    cavityModelCreator creator;

    _cavity = creator.create( bdname );

    

}


/** Destructor */

spotSample::~spotSample() {}



/** Read nodes on boundary */

const void spotSample::readNodes( const string& bdname ) {
    

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

}




/** Read position of nodes on boundary */

const void spotSample::readLocations( const string& bdname ) {

    
    // Look for nodes location

    _location.resize( _nodes.size() );

    for(uint i = 0 ; i < _location.size() ; i++)
	_location[i].resize(3);
    

    ifstream inFile;
    
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



/** Create spots */

const void spotSample::createSpots( const string& bdname ) {


    // Radius model
    
    dictionary dict("properties/adhesiveProperties");

    string rtype = dict.lookUpOrDefault<string>("Boundaries/" + bdname + "/sampleRadiusType", "fixed");
    
    spotRadiusCreator rdcreator;

    std::unique_ptr<spotRadius> rmodel = rdcreator.create(rtype);

    
    
    uint nspots     = (uint)dict.lookUp<scalar>("Boundaries/" + bdname + "/nSpots");

    uint meanRadius = (uint)dict.lookUp<scalar>("Boundaries/" + bdname + "/radius");

    uint devRadius  = (uint)dict.lookUpOrDefault<scalar>("Boundaries/" + bdname + "/dev", 0);

    
    // Assign spots

    _spots.resize(nspots);

    vector<uint> rad = rmodel->radius(nspots, meanRadius, devRadius);

    for(uint i = 0 ; i < nspots ; i++)
	_spots[i] = std::make_pair(0,rad[i]);


}




/** Compute adhesive coefficients */

const std::vector< std::pair<uint, scalar> > spotSample::computeAds() const {

    return _cavity->coeffs(_nodes, _location, _spots);
    
}
