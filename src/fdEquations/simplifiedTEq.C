#include <simplifiedTEq.H>

#include <sphere.H>


using namespace std;


/** Default constructor */

simplifiedTEq::simplifiedTEq( const latticeMesh& mesh, timeOptions& Time, scalarField& T)
    : _mesh(mesh),
      _Time(Time) {



    // Read other constants

    dictionary macroDict("properties/macroProperties");

    _Cv     = macroDict.lookUp<scalar>("Energy/HeatSource/Constants/Cv");

    alpha_1 = macroDict.lookUp<scalar>("Energy/HeatSource/Constants/alpha_1");

    alpha_2 = macroDict.lookUp<scalar>("Energy/HeatSource/Constants/alpha_2");

    _Tau    = macroDict.lookUp< vector<scalar> >( "Energy/LBModel/Tau" );


    

    // Qaux initial values

    for( uint k = 0 ; k < mesh.lmodel()->q() ; k++ )
	Qaux.addElement( 0.5 - 1.0/_Tau[k], k, k );

    Qaux.addElement( (0.5*_Tau[3]-1.0)/_Tau[3], 3, 4 );

    Qaux.addElement( (0.5*_Tau[5]-1.0)/_Tau[5], 5, 6 );    

    

    
    
    // Create eos

    EOSCreator creator;

    eos = creator.create("properties/macroProperties", "Navier-Stokes");
    

    


    // Read boundary conditions for T and set initial values at boundaries

    const map< string, vector<uint> >& bnd = _mesh.boundaries();

    dictionary dict("start/initialFields");

    for( const auto& bd : bnd ) {


	string bdtype = dict.lookUpOrDefault<string>( T.fieldName() + "/boundaryField/" + bd.first + "/type", "none");

	
	if( bdtype == "fixedT" ) {

	    scalar val = dict.lookUp<scalar>( T.fieldName() + "/boundaryField/" + bd.first + "/value");

	    for( const auto& id : bd.second ) {

		T[id] = val;

	    }

	}



	else {

	    if( bdtype == "fixedTSpots" ) {

		scalar val = dict.lookUp<scalar>( T.fieldName() + "/boundaryField/" + bd.first + "/value");

		vector<string> spotList = dict.bracedEntriesNames( T.fieldName() + "/boundaryField/" + bd.first + "/Spots");		


		// Move over spots and assign condition

		for( auto spot : spotList ) {

		

		    // Load type

		    string entry = T.fieldName() + "/boundaryField/" + bd.first + "/Spots/" + spot;

		    const string sptype = dict.lookUp<string>( entry + "/type" );


		    if( sptype == "sphere" ) {

		    	sphere sph( "start/initialFields", entry );

		    	scalar lval = dict.lookUp<scalar>( entry + "/value" );

			for( const auto& id : bd.second ) {

		    	    if(  sph.isInside( _mesh.latticePoint(id) )  ) {

		    		T[id] = lval;

		    	    }

			    else {

				T[id] = val;

			    }

		    	}

		    }

		    else {

		    	cout << " [ERROR] Unrecognized shape type " << sptype << " for fixedTSpots" << endl << endl;

		    	exit(1);

		    }



		}

		

	    }

	}

	

    }
    
}




/** Default destructor */

simplifiedTEq::~simplifiedTEq() {}




/** Source value at specific node */

const scalar simplifiedTEq::heatSource( const uint id ) const {

    return 0;

}




// /** Update scalar field using equation info */

// const void simplifiedTEq::rval( scalarField& field, const scalarField& rho, const vectorField& U, const scalarField& T ) {


//     // Cached variables
    
//     scalar gradRho[3] = {0,0,0};

//     scalar gradT[3] = {0,0,0};

//     scalar lapT(0);
    
    

//     // Move over local points and update field if not on boundary

//     for( uint id = 0 ; id < _mesh.local() ; id++ ) {


// 	// Initialize field value
	
// 	field[id] = 0;
	

// 	if( !_mesh.isOnBoundary(id) ) {


// 	    // \nabla^2 T

// 	    lapT = T.laplacian(id);
	    

	    

//     	    // Convective term: -U \dot \nabla T

// 	    T.grad(gradT, id);

// 	    for( uint j = 0 ; j < 3 ; j++ )
// 		field[id] -= U.at(id)[j] * gradT[j];
		


// 	    // // Only for tests: -T \nabla \cdot U

// 	    // field[id] -= T.at(id) * U.div(id);
	    
	    

// 	    // Diffusive term
	    
// 	    switch( thermal ) {
		
// 	    case thmodel::CD:

		
// 		// Diffusive term: chi \nabla^2 T

// 		field[id] += _lambda * lapT;



// 		// Diffusive term: chi (\nabla rho) \cdot (\nabla T) / \rho
	    
// 		rho.grad(gradRho, id);

// 		for( uint j = 0 ; j < 3 ; j++ )
// 		    field[id] += ( _lambda / rho.at(id) ) * gradRho[j] * gradT[j];

		  		    
	    	   
// 		break;



// 	    case thmodel::CC:

		
// 		// Diffusive term: chi \nabla^2 T

// 		field[id] += _lambda * lapT;
		
// 		break;		


// 	    }




//     	    // Extra term  
	    
// 	    scalar divU( U.div(id) );
	    
//     	    field[id] -= T.at(id) * eos->dp_dT(rho.at(id), T.at(id)) * divU / ( rho.at(id) * _Cv );
	    	   

	    

// 	}



//     }



//     field.sync();
    

// }
