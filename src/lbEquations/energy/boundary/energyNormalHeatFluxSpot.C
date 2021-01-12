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
	    

	    // Create sphere
	    
	    sphere sph( "start/boundaries", entry );	    


	    // Sphere properties
	    
	    const vector<scalar>& centre = sph.centre();

	    const scalar radius = sph.radius();

	    scalar lval  = dict.lookUp<scalar>( entry + "/value" );

	    scalar limit = dict.lookUpOrDefault<scalar>( entry + "/spotLimit", 1000 );

	    scalar eps = dict.lookUpOrDefault<scalar>( entry + "/eps", 1 ); 	    

	    

	    // Assign over nodes
	    
	    for( uint i = 0 ; i < _nodes.size() ; i++ ) {


		uint id = _nodes[i];

		
		// Compute distance to sphere centre

		const vector<int>& point = mesh.latticePoint(id);

		scalar dist = sqrt(  (point[0]-centre[0]) * (point[0]-centre[0])
		                  +  (point[1]-centre[1]) * (point[1]-centre[1])
		                  +  (point[2]-centre[2]) * (point[2]-centre[2]) );



		// Linear interpolation. Piecewise linear
		
		if( dist <= (radius - eps) ) {

		    _grad[i] = lval;

		    _Tlimit[i] = limit;

		}

		else {

		    if( dist <= (radius + eps) ) {

			
			// Temperature gradient distribution
			
			scalar a = ( _grad[i] - lval ) / (2.0*eps);

			scalar b = lval - a * (radius-eps);
			
			_grad[i] = a*dist + b;



			// Temperature limit distribution
			
			a = ( _Tlimit[i] - limit ) / (2.0*eps);

			b = limit - a * (radius-eps);
			
			_Tlimit[i] = a*dist + b;
			

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
