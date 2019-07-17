#include <eqFixedTSpot.H>

using namespace std;


/** Constructor */

eqFixedTSpot::eqFixedTSpot( const std::string& eqName,
				    const std::string& bdName,
				    const latticeMesh& mesh,
				    const scalarField& rho,
				    const scalarField& T,
				    const vectorField& U,
				    pdfField& pdf )

    : eqFixedT( eqName, bdName, mesh, rho, T, U, pdf ) {
    

    // Spots names

    dictionary dict("start/boundaries");
    
    vector<string> spotList = dict.bracedEntriesNames(eqName + "/" + bdName + "/Spots");


    // Move over spots and assign condition

    for( auto spot : spotList ) {
       
	
	// Load type

	string entry = eqName + "/" + bdName + "/Spots/" + spot;

	const string sptype = dict.lookUp<string>( entry + "/type" );


	if( sptype == "sphere" ) {

	    sphere sph( "start/boundaries", entry );

	    scalar lval = dict.lookUp<scalar>( entry + "/value" );

	    for( uint i = 0 ; i < _nodes.size() ; i++ ) {

		uint id = _nodes[i];

		if(  sph.isInside( _mesh.latticePoint(id) )  ) {

		    _bndVal[i] = lval;

		}
		

	    }

	}

	else {

	    cout << " [ERROR] Unrecognized shape type " << sptype << " for fixedTSpots" << endl << endl;

	    exit(1);

	}

    }
    


    
}




/** Destructor */

eqFixedTSpot::~eqFixedTSpot() {}



/** Update pdf field */

void eqFixedTSpot::update( const energyEquation* eeq ) {
    
    eqFixedT::update( eeq );

}
