#include <energyNormalHeatFluxSpot.H>

using namespace std;


/** Constructor */

energyNormalHeatFluxSpot::energyNormalHeatFluxSpot( const std::string& eqName,
						    const std::string& bdName,
						    const latticeMesh& mesh,
						    const scalarField& rho,
						    const scalarField& T,
						    const vectorField& U,
						    pdfField& pdf )

    : energyNormalHeatFlux( eqName, bdName, mesh, rho, T, U, pdf ) {
    

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


	    // Sphere properties
	    
	    const vector<scalar>& centre = sph.centre();

	    const scalar radius = sph.radius();
	    

	    // Assign over nodes
	    
	    for( uint i = 0 ; i < _nodes.size() ; i++ ) {

		uint id = _nodes[i];

		// if(  sph.isInside( _mesh.latticePoint(id) )  ) {

		//     _grad[i] = lval;

		// }

		
		// Compute distance to sphere centre

		const vector<int>& point = mesh.latticePoint(id);

		scalar dist = sqrt(  (point[0]-centre[0]) * (point[0]-centre[0])
		                  +  (point[1]-centre[1]) * (point[1]-centre[1])
		                  +  (point[2]-centre[2]) * (point[2]-centre[2]) );



		// Linear interpolation. Piecewise linear
		
		if( dist <= (radius - 1) ) {

		    _grad[i] = lval;

		}

		else {

		    if( dist <= (radius + 1) ) {

			scalar a = ( _grad[i] - lval ) / 2.0;

			scalar b = lval - a * (radius-1.0);
			
			_grad[i] = a*dist + b;

		    }

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

energyNormalHeatFluxSpot::~energyNormalHeatFluxSpot() {}



/** Update pdf field */

void energyNormalHeatFluxSpot::update( const energyEquation* eeq ) {
    
    energyNormalHeatFlux::update( eeq );

}
