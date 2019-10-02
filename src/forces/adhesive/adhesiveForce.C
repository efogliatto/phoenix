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

	    if( bd.second.size() != 0 ) {

		uint nid;

		scalar g;
	    
		for( const auto &id : bd.second ) {

		    infile >> nid;

		    infile >> g;

		    _Gads[id] = g;

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
