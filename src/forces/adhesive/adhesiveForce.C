#include <adhesiveForce.H>

using namespace std;


/** Constructor */

adhesiveForce::adhesiveForce( const string& dictName,
			      const string& eqName,
			      const latticeMesh& mesh,
			      const interactionForce* Fi,
			      timeOptions& Time)

    : _mesh(mesh),
      _Fi(Fi),
      _force(mesh, Time, "Fads", IO::NO_READ, IO::NO_WRITE) {
    

    // Initialize force

    for( uint i = 0 ; i < _mesh.local() ; i++ ) {

	_force[i][0] = 0;
	_force[i][1] = 0;
	_force[i][2] = 0;	

    }
    
    

    // Read adsorption coefficient for each node on boundary

    const map< string, vector<uint> >& boundary = _mesh.boundaries();

        
    for( const auto &bd : boundary ) {



	// If file is present, read from file

	ifstream infile( "processor" + to_string(_mesh.pid()) + "/lattice/" + bd.first + "_gads" );

	if( infile.good() ) {

	    uint npoints(0);

	    scalar g;

	    infile >> npoints;
	    
	    
	    if( npoints > 0 ) {

		uint nid(0);
	    
	    	for( uint k = 0 ; k < npoints ; k++ ) {

	    	    infile >> nid;

	    	    infile >> g;

	    	    _Gads[nid] = g;

	    	}

	    }

	}

	else {
	    

	    dictionary dict("start/boundaries");
		
	    scalar g_ads = dict.lookUpOrDefault<scalar>( eqName + "/" + bd.first + "/Gads", 0 );

	    if(  ( g_ads != 0 )  &&  ( bd.second.size() != 0 )  ) {
	    
		for( const auto &id : bd.second ) {

		    _Gads[id] = g_ads;

		}

	    }


	}

    }
    
    
    
}


/** Destructor */

adhesiveForce::~adhesiveForce() {}
