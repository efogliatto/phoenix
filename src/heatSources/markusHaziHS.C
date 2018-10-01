#include <markusHaziHS.H>

using namespace std;


/** Constructor */

markusHaziHS::markusHaziHS( const string& dictName,
			    const string& eqName,
			    const latticeMesh& mesh,
			    timeOptions& Time )

    : heatSource(dictName, eqName, mesh, Time) {


    
    // Create eos

    EOSCreator creator;

    eos = creator.create(dictName, "Navier-Stokes");

    
    // Read coefficients from dictionary

    dictionary dict( dictName );

    _Tau = dict.lookUp< vector<scalar> >( eqName + "/LBModel/Tau" );

    _Cv = dict.lookUp<scalar>( eqName + "/HeatSource/Constants/Cv" );

    _a1 = dict.lookUp<scalar>( eqName + "/HeatSource/Constants/alpha_1" );

    _a2 = dict.lookUp<scalar>( eqName + "/HeatSource/Constants/alpha_2" );
    

}




/** Update source field */

void markusHaziHS::update( const scalarField& rho, const scalarField& T, const vectorField& U ) {


    // Constants

    scalar gradT[3]   = {0,0,0};

    scalar gradRho[3] = {0,0,0};

    const uint np = _mesh.local();



    // Move over points

    for( uint id = 0 ; id < np ; id++ ) {


    
	// Thermal difusivity
    
	scalar lambda;

	if( _mesh.lmodel()->name() == "D2Q9" ) {

	    lambda = rho.at(id) * _Cv * (1/_Tau[3] - 0.5) * (4.0 + 3.0 * _a1  + 2.0 * _a2) / 6.0;	

	}

	else {

	    cout << " [ERROR]  Heat source model not implemented" << endl << endl;

	    exit(1);

	}


	// Cached scalar values

	const scalar _rho = rho.at(id);

	const scalar _T = T.at(id);



	// Scalar fields gradients

	T.grad(gradT, id);

	rho.grad(gradRho, id);

	
	scalar first(0);

	for(uint j = 0 ; j < 3 ; j++)
	    first += gradT[j] * gradRho[j];

	first = first * lambda / (_rho * _rho * _Cv);


	
	// Velocity divergence term

	scalar second = U.div(id) * _T * ( 1.0   -   eos->dp_dT(_rho, _T) / (_rho * _Cv) );


	
	// Update source at node
	
	_source[id] = first + second;
       
    
    }



    // Sync source

    _source.sync();
    

}
