#include <rndSpots.H>

using namespace std;


/** Constructor */

rndSpots::rndSpots( const string& dictName,
		    const string& eqName,
		    const latticeMesh& mesh,
		    const interactionForce* Fi,
		    timeOptions& Time )

    : phiBasedMod(dictName, eqName, mesh, Fi, Time) {


    
    // Update adhesive constant distribution

    const map< string, vector<uint> >& boundary = _mesh.boundaries();

    dictionary dict("start/boundaries");
        
    for( const auto &bd : boundary ) {

	int nspots = dict.lookUpOrDefault<int>( eqName + "/" + bd.first + "/nspots", 0 );

	
	
    	// scalar g_ads = dict.lookUpOrDefault<scalar>( eqName + "/" + bd.first + "/Gads", 0 );

    	// if(  ( g_ads != 0 )  &&  ( bd.second.size() != 0 )  ) {
	    
	//     for( const auto &id : bd.second ) {

	// 	_Gads[id] = g_ads;

	//     }

	// }

    }    

}



/** Destructor */

rndSpots::~rndSpots() {}



/** Update force field */

void rndSpots::update( scalarField& rho, scalarField& T ) {

    phiBasedMod::update(rho, T);

}
