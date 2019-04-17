#include <TEquation.H>

#include <sphere.H>


using namespace std;


/** Default constructor */

TEquation::TEquation( const latticeMesh& mesh, timeOptions& Time, scalarField& T)
    : _mesh(mesh),
      _Time(Time) {



    // Read other constants

    dictionary macroDict("properties/macroProperties");

    _Cv = macroDict.lookUp<scalar>("Energy/Cv");

    _lambda = macroDict.lookUp<scalar>("Energy/lambda");    


    string dmodel = macroDict.lookUp<string>("Energy/ConductivityModel");

    if( dmodel == "constDiff" ) {

	thermal = thmodel::CD;

    }

    else {

	if( dmodel == "constCond" ) {

	    thermal = thmodel::CC;

	}

    }



    
    // Create eos

    EOSCreator creator;

    eos = creator.create("properties/macroProperties", "Navier-Stokes");
    

    


    // Read boundary conditions for T and set initial values at boundaries

    const map< string, vector<uint> >& bnd = _mesh.boundaries();

    dictionary dict("start/initialFields");

    for( const auto& bd : bnd ) {


	string bdtype = dict.lookUpOrDefault<string>("T/boundaryField/" + bd.first + "/type", "none");

	
	if( bdtype == "fixedT" ) {

	    scalar val = dict.lookUp<scalar>("T/boundaryField/" + bd.first + "/value");

	    for( const auto& id : bd.second ) {

		T[id] = val;

	    }

	}



	else {

	    if( bdtype == "fixedTSpots" ) {

		scalar val = dict.lookUp<scalar>("T/boundaryField/" + bd.first + "/value");

		vector<string> spotList = dict.bracedEntriesNames("T/boundaryField/" + bd.first + "/Spots");		


		// Move over spots and assign condition

		for( auto spot : spotList ) {

		

		    // Load type

		    string entry = "T/boundaryField/" + bd.first + "/Spots/" + spot;

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

TEquation::~TEquation() {}




/** Update scalar field using equation info */

const void TEquation::rval( scalarField& field, const scalarField& rho, const vectorField& U, const scalarField& T ) {


    // Lattice constants
    
    const uint q( _mesh.lmodel()->q() );

    const vector< vector<int> >& nb = _mesh.nbArray();

    const vector<uint>& reverse = _mesh.lmodel()->reverse();

    const vector< vector<int> >& vel = _mesh.lmodel()->lvel();

    const scalar cs2 = _mesh.lmodel()->cs2();

    const vector<scalar>& omega = _mesh.lmodel()->omega();

    scalar gradRho[3] = {0,0,0};
    
    

    // Move over local points and update field if not on boundary

    for( uint id = 0 ; id < _mesh.local() ; id++ ) {


	// Initialize field value
	
	field[id] = 0;
	

	if( !_mesh.isOnBoundary(id) ) {
	

    	    // Convective term: -U \dot \nabla T

    	    for( uint k = 1 ; k < q ; k++ ) {

    	    	int nbid = nb[id][reverse[k]];

		for( uint j = 0 ; j < 3 ; j++ )
		    field[id] -= T.at(nbid) * omega[k] * vel[k][j] * U.at(id)[j] / cs2;		

    	    }


	    // Diffusive term
	    
	    switch( thermal ) {
		
	    case thmodel::CD:

		
		// Diffusive term: chi \nabla^2 T

		for( uint k = 1 ; k < q ; k++ ) {

		    int nbid = nb[id][reverse[k]];
		
		    field[id] += 2 * omega[k] * _lambda * ( T.at(nbid) - T.at(id) ) / cs2;

		}



		// Diffusive term: chi (\nabla rho) \cdot (\nabla T) / \rho
	    
		rho.grad(gradRho, id);
	    	   
		for( uint k = 1 ; k < q ; k++ ) {

		    int nbid = nb[id][reverse[k]];

		    for( uint j = 0 ; j < 3 ; j++ )    	    	   
			field[id] += _lambda * omega[k] * vel[k][j] * gradRho[j] * T.at(nbid) / (cs2 * rho.at(id));

		}

	    


		break;



	    case thmodel::CC:

		
		// Diffusive term: chi \nabla^2 T

		for( uint k = 1 ; k < q ; k++ ) {

		    int nbid = nb[id][reverse[k]];
		
		    field[id] += 2 * omega[k] * _lambda * ( T.at(nbid) - T.at(id) ) / (rho.at(id) * cs2 * _Cv);

		}   


		break;		


	    }






    	    // Extra term

    	    // scalar divU = 0.5 * U.at( nb[id][3] )[0]
    	    // 	        - 0.5 * U.at( nb[id][1] )[0]
    	    // 	        + 0.5 * U.at( nb[id][4] )[1]
    	    // 	        - 0.5 * U.at( nb[id][2] )[1];	  
	    
	    scalar divU( U.div(id) );
	    
    	    field[id] -= T.at(id) * eos->dp_dT(rho.at(id), T.at(id)) * divU / ( rho.at(id) * _Cv );
	    	   

	    

	}



    }



    field.sync();
    

}
