#include <GongSRTEq.H>

#include <algebra.H>

using namespace std;


/** Default constructor */

GongSRTEq::GongSRTEq( const std::string& name,
		      const latticeMesh& mesh_,
		      timeOptions& Time_,
		      pdfField& pdf_,
		      const scalarField& rho_,
		      const vectorField& U_,
		      scalarField& T_)

    : energyEquation(name, mesh_, Time_, pdf_, rho_, U_, T_) {


    // Read constants

    dictionary dict("properties/macroProperties");

    _kappa = dict.lookUp<scalar>( name + "/HeatSource/Constants/kappa");

    _Cv = dict.lookUp<scalar>( name + "/HeatSource/Constants/Cv");
    

    
    // Create eos

    EOSCreator creator;

    eos = creator.create("properties/macroProperties", "Navier-Stokes");

    

    // Set equilibrium pdf values only if time is start time

    if( Time.currentTime() == Time.timeList()[0] )  {

	GongSRTEq::setEquilibrium();

    }


    // Read conductivity model
    
    string th = dict.lookUp<string>( name + "/HeatSource/ConductivityModel" );

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



/** Default destructor */

GongSRTEq::~GongSRTEq() {}




/** Equilibrium in population space */

const void GongSRTEq::eqPS( vector<scalar>& n, const uint& id ) const {

    GongSRTEq::eqPS( n, T.at(id), U.at(id), 0 );
    
}



/** Equilibrium in population space */

const void GongSRTEq::eqPS( vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_, const scalar& hs ) const {


    // Lattice constants
    
    const uint q = mesh.lmodel()->q();

    const scalar cs2 = mesh.lmodel()->cs2();

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();

    const vector<scalar>& omega = mesh.lmodel()->omega();



    // Move over velocities
    
    for( uint k = 0 ; k < q ; k++ ) {

    	scalar alpha(0),
    	    beta(0);

	
    	// Dot product
	
    	for( uint j = 0 ; j < 3 ; j++ ) {

    	    alpha += vel[k][j] * U_[j];

    	    beta += U_[j] * U_[j];

    	}

    	n[k] = T_ * omega[k] * ( 1 + alpha/cs2   +   0.5 * alpha * alpha / (cs2*cs2)  -  0.5 * beta / cs2 );
	
    }
    
    
    
}



/** Set pdf to equilibrium values */

const void GongSRTEq::setEquilibrium() {

    const scalar q = mesh.lmodel()->q();
    
    vector<scalar> n( q );
    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

    	eqPS(n,i);

    	for(uint k = 0 ; k < q ; k++) 
    	    _pdf.set(i,k,n[k]);	    

    }

}



/** Collision process */

const void GongSRTEq::collision() {


    // Lattice constants
    
    const uint nodes = mesh.local();

    const uint q = mesh.lmodel()->q();

    const vector<scalar>& omega = mesh.lmodel()->omega();


    // Partial distributions
    
    vector<scalar> g_eq(q);
    
    
    
    // Move over local points

    for( uint id = 0 ; id < nodes ; id++ ) {


	// Compute equilibrium
	
	GongSRTEq::eqPS(g_eq, id);


	// Local heat source

	scalar heat(0);

	{

	    scalar _rho = rho.at(id);

	    scalar first(0);

	    scalar gradT[3]   = {0,0,0};

	    scalar gradRho[3] = {0,0,0};



	    // Temperature laplacian

	    scalar lapT = T.laplacian(id);
	


	    // Compute first term 

	    switch( thermal ) {

	    case thmodel::constCond:

		first = lapT * ( _kappa / ( _rho *_Cv) - (_Tau[0]-0.5) / 3.0 );
	    
		break;


	    case thmodel::constDiff:
	    
		first = 0;
	    
		T.grad(gradT, id);

		rho.grad(gradRho, id);



		// Dot product
	    
		for(uint j = 0 ; j < 3 ; j++)
		    first += gradT[j] * gradRho[j];

		first = (_kappa / _rho) * ( first +  _rho * lapT )   -   (_Tau[0]-0.5) * lapT / 3.0;
	    
	    
		break;

	    }
	
	
	    // Velocity divergence term

	    scalar second = U.div(id) * T.at(id) * ( 1.0   -   eos->dp_dT(_rho, T.at(id)) / (_rho * _Cv) );


	    
	    heat = first + second;

	}


	
	// Collide in population space

	for( uint k = 0 ; k < q ; k++ ) {

	    _pdf[id][k] = _pdf[id][k]
		        - (1.0/_Tau[0]) * (_pdf[id][k] - g_eq[k])
		        + omega[k] * heat;

	}

	

    }


    // Sync ghost nodes
    
    _pdf.sync();

}







/** Update macroscopic temperature */

const void GongSRTEq::updateMacroTemperature() {
    
    for( uint i = 0 ; i < mesh.npoints() ; i++ )
    	T[i] = energyEquation::localTemperature(i);

}




/** Thermal conductivity at node */

const scalar GongSRTEq::thermalCond( const uint& id ) const {

    scalar k(0);

    return k;

}



/** Diffusivity constant recovered at macroscopic level */

const scalar GongSRTEq::diffusivityConstant() const {

    return (_Tau[0] - 0.5) / 3.0;

}
