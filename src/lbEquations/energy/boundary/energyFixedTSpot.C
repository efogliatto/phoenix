#include <energyFixedTSpot.H>

using namespace std;


/** Constructor */

energyFixedTSpot::energyFixedTSpot( const std::string& eqName,
				    const std::string& bdName,
				    const latticeMesh& mesh,
				    const scalarField& rho,
				    const scalarField& T,
				    const vectorField& U,
				    pdfField& pdf )

    : energyFixedT( eqName, bdName, mesh, rho, T, U, pdf ) {
    

    // Spots names

    dictionary dict("start/boundaries");
    
    vector<string> spotList = dict.bracedEntriesNames(eqName + "/" + bdName + "/Spots");


    // Move over spots and assign condition

    for( auto spot : spotList ) {

	
	// Load type

	const string sptype = dict.lookUp<string>( eqName + "/" + bdName + "/Spots/" + spot + "/type" );


	if( sptype == "sphere" ) {

	    sphere sph( "start/boundaries", eqName + "/" + bdName + "/Spots/" + spot );

	    

	    for( uint i = 0 ; i < _nodes.size() ; i++ ) {

		uint id = _nodes[i];	    
		

	    }

	}

	else {

	    cout << " [ERROR] Unrecognized shape type" << endl << endl;

	    exit(1);

	}

    }
    


    
}




/** Destructor */

energyFixedTSpot::~energyFixedTSpot() {}



/** Update pdf field */

void energyFixedTSpot::update( const energyEquation* eeq ) {
    
    energyFixedT::update( eeq );

}
