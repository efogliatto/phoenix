#include <energyFixedT.H>

#include <random>

using namespace std;



/** Constructor */

energyFixedT::energyFixedT( const std::string& eqName,
			    const std::string& bdName,
			    const latticeMesh& mesh,	  
			    const scalarField& rho,
			    const scalarField& T,
			    const vectorField& U,
			    pdfField& pdf )
    
    : energyBndCond(mesh, rho, T, U, pdf, bdName) {

    

    // Load boundary value

    dictionary dict("start/boundaries");

    scalar val = dict.lookUp<scalar>( eqName + "/" + bdName + "/value" );


    // Resize values at boundary and assign

    _bndVal.resize( _nodes.size() );

    for( uint i = 0 ; i < _bndVal.size() ; i++ )
	_bndVal[i] = val;



    // Perturbation

    _pert = dict.lookUpOrDefault<scalar>( eqName + "/" + bdName + "/perturbation", 0 );





    // Resize neighbour indices, compute normals and check indexing

    _nbid.resize( _nodes.size() );

    const vector< vector<int> >& nb = mesh.nbArray();

    const vector<uint>& reverse = mesh.lmodel()->reverse();

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();

    const uint q = mesh.lmodel()->q();


    
    for( uint i = 0 ; i < _nodes.size() ; i++ ) {

	uint nid = _nodes[i];


	// Check which is the normal direction, using unknown neighbours

	scalar normal[3] = {0,0,0};

	for( uint k = 0 ; k < q ; k++ ) {

	    if( nb[nid][ reverse[k] ] == -1 ) {

		for( uint j = 0 ; j < 3 ; j++ ) {

		    normal[j] += (scalar)vel[k][j];

		}

	    }

	}

	scalar nmag(0);

	for( uint j = 0 ; j < 3 ; j++ )
	    nmag += normal[j] * normal[j];
	    


	if( nmag!=0 ) {

	    nmag = sqrt(nmag);

	    for( uint j = 0 ; j < 3 ; j++ ) {

		normal[j] = normal[j] / nmag;

	    }

	    _normal.push_back( {(int)normal[0], (int)normal[1], (int)normal[2]} );
	    
	}

	else {

	    cout << " [ERROR]  Unable to detect normal on node " << i << " over boundary " << bdName;

	    exit(1);

	}
             

	// Detect correspondence with lattice velocities

	int ln = -1;

	for( uint k = 0 ; k < q ; k++ ) {
	    
	    if(      ( vel[k][0] == (int)normal[0] )
		 &&  ( vel[k][1] == (int)normal[1] )
		 &&  ( vel[k][2] == (int)normal[2] )  ) {
		
	    	ln = k;

	    }

	}

	if(ln == -1) {

	    cout << " [ERROR]  Unable to detect normal on node " << _nodes[i] << " over boundary " << bdName;

	    exit(1);

	}

	else {

	    _nbid[i] = nb[nid][ln];

	}
	

    }    
    

}



/** Destructor */

energyFixedT::~energyFixedT() {}



/** Update pdf field. NEBB */

void energyFixedT::update( const energyEquation* eeq ) {



    // // Lattice constants
    
    // const uint q = _mesh.lmodel()->q();

    // const vector< vector<int> >& nb = _mesh.nbArray();

    // const vector<uint>& reverse = _mesh.lmodel()->reverse();     

    // vector<scalar> f_eq_bnd(q);

    // vector<scalar> Uw = {0,0,0};


    

    // // Random seeds
    
    // uniform_real_distribution<scalar> unif( (100.0-_pert)/100.0, (100.0+_pert)/100.0);

    // default_random_engine re;

    // re.seed( _mesh.pid() );

    


    // // Move over boundary elements

    // for( uint i = 0 ; i < _nodes.size() ; i++ ) {
	

    // 	uint id = _nodes[i];

    // 	scalar Tw = unif(re) * _bndVal[i];



    // 	// Density and velocity at neighbour node

    // 	for(uint j = 0 ; j < 3 ; j++)			
    // 	    Uw[j] = _U.at(id,j);


	

    // 	// Equilibrium populations over boundary

    // 	eeq->eqPS( f_eq_bnd, Tw, Uw, 0 );


	
    // 	// Update unknowk distributions

    // 	for( uint k = 1 ; k < q ; k++ ) {	    	    		       		   		    
			
    // 	    if( nb[id][k] == -1 ) {

    // 		// _pdf[id][k] = _pdf[id][reverse[k]];
    // 		_pdf[id][k] = f_eq_bnd[k];

    // 	    }

    // 	}


    // 	// Beta constant

    // 	scalar beta(0);

    // 	scalar kn(0);

    // 	scalar unk(0);


    // 	for( uint k = 0 ; k < q ; k++ ) {	    	    		       		   		    
			
    // 	    if( nb[id][k] == -1 ) {

    // 		unk += _pdf[id][k];

    // 	    }

    // 	    else {

    // 		kn += _pdf[id][k];

    // 	    }

    // 	}

    // 	beta = (Tw - kn) / unk;



    // 	for( uint k = 1 ; k < q ; k++ ) {	    	    		       		   		    
			
    // 	    if( nb[id][k] == -1 ) {
		
    // 		_pdf[id][k] = beta * _pdf[id][k];

    // 	    }

    // 	}
	
	

    // }

    

}





// /** Update pdf field */

// // Non-equilibrium extrapolation

// void energyFixedT::update( const energyEquation* eeq ) {



//     // Lattice constants
    
//     const uint q = _mesh.lmodel()->q();

//     vector<scalar> f_eq_nb(q);

//     vector<scalar> f_eq_bnd(q);

//     vector<scalar> Unbid = {0,0,0};

//     vector<scalar> Uw = {0,0,0};


    

//     // Random seeds
    
//     uniform_real_distribution<scalar> unif( (100.0-_pert)/100.0, (100.0+_pert)/100.0);

//     default_random_engine re;

//     re.seed( _mesh.pid() );

    


//     // Move over boundary elements

//     for( uint i = 0 ; i < _nodes.size() ; i++ ) {
	

//     	uint id = _nodes[i];

// 	scalar Tw = unif(re) * _bndVal[i];


	
// 	// Density and velocity at neighbour node

// 	const uint nbid = _nbid[i];
	


// 	// Velocity at boundary

// 	for(uint j = 0 ; j < 3 ; j++) {
			
// 	    Unbid[j] = _U.at(nbid,j);

// 	    // Uw[j] = _U.at(id,j);

// 	}

	

// 	// Equilibrium

// 	eeq->eqPS( f_eq_nb, _T.at(nbid), Unbid, 0 );

// 	eeq->eqPS( f_eq_bnd, Tw, Uw, 0 );

// 	// eeq->eqPS( f_eq_nb, _T.at(nbid), Unbid, eeq->heat(id) );

// 	// eeq->eqPS( f_eq_bnd, Tw, Uw, eeq->heat(id) );
	


	
// 	// Update distribution

//     	for( uint k = 0 ; k < q ; k++ ) {	    	    		       		   		    
			
// 	    _pdf[id][k] = f_eq_bnd[k] + _pdf[nbid][k] - f_eq_nb[k];

//     	}
	

//     }

    

// }







// /** Update pdf field */

// // Equiilbrium scheme

// void energyFixedT::update( const energyEquation* eeq ) {



//     // Lattice constants
    
//     const uint q = _mesh.lmodel()->q();

//     vector<scalar> f_eq_bnd(q);

//     vector<scalar> Uw = {0,0,0};




//     // Move over boundary elements

//     for( uint i = 0 ; i < _nodes.size() ; i++ ) {
	     	

// 	// Equilibrium

//     	uint id = _nodes[i];

// 	scalar Tw = _bndVal[i];
	
// 	eeq->eqPS( f_eq_bnd, Tw, Uw, 0 );	


	
// 	// Update distribution

//     	for( uint k = 0 ; k < q ; k++ ) {	    	    		       		   		    
			
// 	    _pdf.set(id, k, f_eq_bnd[k] );

//     	}
	

//     }

    

// }
