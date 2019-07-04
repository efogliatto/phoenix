#include <simplifiedTEq.H>

#include <sphere.H>


using namespace std;


/** Default constructor */

simplifiedTEq::simplifiedTEq( const latticeMesh& mesh, timeOptions& Time, scalarField& T)
    : _mesh(mesh),
      _Time(Time),
      Tnew( mesh, Time, "Tnew", IO::NO_READ, IO::NO_WRITE ) {



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




/** Update scalar field using equation info */

const void simplifiedTEq::predictor( scalarField& Tstar, const scalarField& rho, const vectorField& U, const scalarField& T ) {

    
    // Lattice constants

    const vector< vector<int> >& nb = _mesh.nbArray();
	    
    const scalarMatrix& invM = _mesh.lmodel()->MRTInvMatrix();

    const uint q = _mesh.lmodel()->q();

    vector<scalar> n_eq(q);

    vector<scalar> n(q);    



    // Move over local points only

    for( uint id = 0 ; id < _mesh.local() ; id++ ) {


	if( !_mesh.isOnBoundary(id) ) {


    	    for( uint k = 0 ; k < q ; k++ ) {

    		int nbid = nb[id][k];

    		if(nbid != -1)
    		    n_eq[k] = simplifiedTEq::neq( T, U, nbid, k );

    	    }


    	    invM.matDotVec(n_eq, n);


    	    Tstar[id] = 0.5*simplifiedTEq::heatSource(rho, T, U, id);

    	    for( uint k = 0 ; k < q ; k++ )
    		Tstar[id] += n[k];

    	}

	else {

	    Tstar[id] = T.at(id);

	}

		
    }


    Tstar.sync();
    

}





/** Corrector step (diffusion) */

const void simplifiedTEq::corrector( scalarField& T, const scalarField& rho, const vectorField& U, const scalarField& Tstar ) {


    // // Lattice constants

    // const vector< vector<int> >& nb = _mesh.nbArray();
	    
    // const scalarMatrix& invM = _mesh.lmodel()->MRTInvMatrix();

    // vector<uint> reverse = _mesh.lmodel()->reverse();	    

    // const uint q = _mesh.lmodel()->q();

    // vector<scalar> n_eq(q);

    // vector<scalar> n(q);
    
    
    
    // for( uint id = 0 ; id < _mesh.local() ; id++ ) {

	
    // 	if( !_mesh.isOnBoundary(id) ) {
	    

    // 	    for( uint k = 0 ; k < q ; k++ ) {

    // 		int nbid = nb[id][k];

    // 		int nbplus = nb[id][reverse[k]];

    // 		if(nbid != -1) {

    // 		    if (nbplus != -1) {
			
    // 			n_eq[k] = simplifiedTEq::neq( Tstar, U, nbplus, k )
    // 			        - simplifiedTEq::neq( Tstar, U, id, k )
    // 			        + simplifiedTEq::neq( T, U, nbid, k )
    // 			        - simplifiedTEq::neq( T, U, id, k );

    // 		    }
		    
    // 		}


    // 	    }


    // 	    Qaux.matDotVec(n_eq, n);

    // 	    invM.matDotVec(n, n_eq);

		   

    // 	    // Update new temperature field

    // 	    Tnew[id] = Tstar.at(id);

    // 	    for( uint k = 0 ; k < q ; k++ )
    // 	    	Tnew[id] += n_eq[k];

    // 	}
		
    // }



    // // Update T values

    // for( uint id = 0 ; id < _mesh.local() ; id++ ) {
	
    // 	if( !_mesh.isOnBoundary(id) )
    // 	    T[id] = Tnew.at(id);

    // }

    // T.sync();









    
    // Update T values

    scalar chi(0);

    switch( _mesh.lmodel()->q() ) {

    case 9:

    	chi = (1/_Tau[3] - 0.5) * (4.0 + 3.0 * alpha_1  + 2.0 * alpha_2) / 6.0;

    	break;


    case 15:

    	chi = (1/_Tau[3] - 0.5) * (6.0 + 11.0 * alpha_1  + alpha_2) / 9.0;

    	break;

    }



    for( uint id = 0 ; id < _mesh.local() ; id++ ) {

	if( !_mesh.isOnBoundary(id) )
	    Tnew[id] = 0.5 * chi * T.laplacian(id);

    }


    for( uint id = 0 ; id < _mesh.local() ; id++ ) {

	if( !_mesh.isOnBoundary(id) )
	    T[id] = Tstar.at(id) + Tnew[id];

    }

        

    T.sync();
    
    
}









/** Source value at specific node (real) */ 

const scalar simplifiedTEq::heatSource( const scalarField& rho, const scalarField& T, const vectorField& U, const uint id ) const {

    
    // Constants

    scalar gradT[3]   = {0,0,0};

    scalar gradRho[3] = {0,0,0};



    // Cached scalar values

    const scalar _rho = rho.at(id);

    const scalar _T = T.at(id);
	
	
    
    // Thermal difusivity
    
    scalar chi(0);

    switch( _mesh.lmodel()->q() ) {

    case 9:

    	chi = (1/_Tau[3] - 0.5) * (4.0 + 3.0 * alpha_1  + 2.0 * alpha_2) / 6.0;

    	break;


    case 15:

    	chi = (1/_Tau[3] - 0.5) * (6.0 + 11.0 * alpha_1  + alpha_2) / 9.0;

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


	
    return first + second;


}





/** Equilibrium value at specific node and velocity */

const scalar simplifiedTEq::neq( const scalarField& T, const vectorField& U, const uint id, const uint vid ) {

    
    const uint q = _mesh.lmodel()->q();

    scalar n(0);
    

    switch(q) {
	

    case 9:
	
	switch(vid) {

	case 0:

	    n = T.at(id);

	    break;


	case 1:

	    n = alpha_1 * T.at(id);

	    break;


	case 2:

	    n = alpha_2 * T.at(id);

	    break;


	case 3:

	    n = T.at(id) * U.at(id)[0];

	    break;


	case 4:

	    n = -T.at(id) * U.at(id)[0];

	    break;


	case 5:

	    n = T.at(id) * U.at(id)[1];

	    break;


	case 6:

	    n = -T.at(id) * U.at(id)[1];

	    break;


	case 7:

	    n = 0;

	    break;


	case 8:

	    n = 0;

	    break;	    

	}
	
	
	break;



    default:

    	cout << " [ERROR]  Equilibrium model not implemented" << endl;

    	exit(1);

	break;	

    }


    return n;

}
