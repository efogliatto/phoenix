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


    // Lattice type

    latticeModel::latticeType ltype = _mesh.lmodel()->type();
    

    // Move over points

    for( uint id = 0 ; id < np ; id++ ) {


	// Cached scalar values

	const scalar _rho = rho.at(id);

	const scalar _T = T.at(id);
	
	
    
	// Thermal difusivity
    
	scalar chi(0);

	switch( ltype ) {


	case latticeModel::latticeType::D2Q9:

	    chi = (1/_Tau[3] - 0.5) * (4.0 + 3.0 * _a1  + 2.0 * _a2) / 6.0;

	    break;


	case latticeModel::latticeType::D3Q15:

	    chi = (1/_Tau[3] - 0.5) * (6.0 + 11.0 * _a1  + _a2) / 9.0;

	    break;

	}

	      



	// Scalar fields gradients

	T.grad(gradT, id);

	rho.grad(gradRho, id);

	
	scalar first(0);

	for(uint j = 0 ; j < 3 ; j++)
	    first += gradT[j] * gradRho[j];

	first = first * chi / _rho;

	
	
	// Velocity divergence term

	scalar second = U.div(id) * _T * ( 1.0   -   eos->dp_dT(_rho, _T) / (_rho * _Cv) );


	
	// Update source at node
	
	_source[id] = first + second;
       
    
    }



    // Sync source

    _source.sync();
    

}




/** Update source field with external temperature gradient */
    
void markusHaziHS::update( const scalarField& rho, const scalarField& T, const vectorField& U, const vectorField& Tgrad ) {


    // Constants

    scalar gradRho[3] = {0,0,0};

    const uint np = _mesh.local();


    // Lattice type

    latticeModel::latticeType ltype = _mesh.lmodel()->type();
    

    // Move over points

    for( uint id = 0 ; id < np ; id++ ) {


	// Cached scalar values

	const scalar _rho = rho.at(id);

	const scalar _T = T.at(id);
	
	
    
	// Thermal difusivity
    
	scalar chi(0);

	switch( ltype ) {


	case latticeModel::latticeType::D2Q9:

	    chi = (1/_Tau[3] - 0.5) * (4.0 + 3.0 * _a1  + 2.0 * _a2) / 6.0;

	    break;


	case latticeModel::latticeType::D3Q15:

	    chi = (1/_Tau[3] - 0.5) * (6.0 + 11.0 * _a1  + _a2) / 9.0;

	    break;

	}

	      



	// Scalar fields gradients

	rho.grad(gradRho, id);

	
	scalar first(0);

	for(uint j = 0 ; j < 3 ; j++)
	    first += Tgrad.at(id)[j] * gradRho[j];

	first = first * chi / _rho;

	
	
	// Velocity divergence term

	scalar second = U.div(id) * _T * ( 1.0   -   eos->dp_dT(_rho, _T) / (_rho * _Cv) );


	
	// Update source at node
	
	_source[id] = first + second;
       
    
    }



    // Sync source

    _source.sync();
    

}
