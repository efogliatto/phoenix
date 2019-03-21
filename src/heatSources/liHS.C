#include <liHS.H>

using namespace std;


/** Constructor */

liHS::liHS( const string& dictName,
	    const string& eqName,
	    const latticeMesh& mesh,
	    timeOptions& Time )

    : heatSource(dictName, eqName, mesh, Time),
      _Lambda( mesh, Time, "lambda", IO::NO_READ, IO::NO_WRITE ) {


    
    // Create eos

    EOSCreator creator;

    eos = creator.create(dictName, "Navier-Stokes");

    
    // Read coefficients from dictionary

    dictionary dict( dictName );

    _Tau = dict.lookUp< vector<scalar> >( eqName + "/LBModel/Tau" );

    _Cv = dict.lookUp<scalar>( eqName + "/HeatSource/Constants/Cv" );

    _kappa = dict.lookUp<scalar>( eqName + "/HeatSource/Constants/kappa" );


    
    // Read conductivity model
    
    string th = dict.lookUp<string>( eqName + "/HeatSource/ConductivityModel" );

    if(th == "constCond") {

	thermal = thmodel::constCond;

    }

    else {

	if(th == "constDiff") {

	    thermal = thmodel::constDiff;

	}

	else {

	    cout << "Undefined thermal conductivity model" << endl;

	    exit(1);

	}

    }
    

}




/** Update source field */

void liHS::update( const scalarField& rho, const scalarField& T, const vectorField& U ) {


    // Constants

    const uint np = _mesh.local();

    const scalar _k = ( 1/_Tau[3] - 0.5 ) * _mesh.lmodel()->cs2();

    

    // Update thermal conductivity field

    if ( thermal == thmodel::constDiff ) {

	for( uint id = 0 ; id < _mesh.npoints() ; id++ ) {

	    _Lambda[id] = _kappa * rho.at(id) * _Cv;

	}

    }

    

    // Move over points

    for( uint id = 0 ; id < np ; id++ ) {


    	// Cached scalar values

    	const scalar _rho = rho.at(id);

    	const scalar _T = T.at(id);

	scalar first(0);

	scalar gradT[3]   = {0,0,0};

	scalar gradLambda[3] = {0,0,0};	



	// Temperature laplacian

	scalar lapT = T.laplacian(id);
	


	// Compute first term 
	
	switch( thermal ) {

	case thmodel::constCond:

	    first = lapT * ( _kappa / ( _rho *_Cv) - _k );
	    
	    break;


	case thmodel::constDiff:

	    T.grad(gradT, id);

	    _Lambda.grad(gradLambda, id);


	    // Dot product
	    
	    for(uint j = 0 ; j < 3 ; j++)
		first += gradT[j] * gradLambda[j];


	    first = ( first + _Lambda[id] * lapT ) / ( _rho * _Cv )   -   _k*lapT;

	    
	    break;

	}
	
	
	
	
    	// Velocity divergence term

    	scalar second = U.div(id) * _T * ( 1.0   -   eos->dp_dT(_rho, _T) / (_rho * _Cv) );



	
    	// Update source at node
	
    	_source[id] = 1.5 * (first + second) - 0.5 * _source[id];


	
	// _source[id] = first + second;
	
    	// _source[id] = 1.5 * (first + second) - 0.5 * _source[id];
       
    
    }



    // Sync source

    _source.sync();
    

}
