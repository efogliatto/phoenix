#include <ppFixedU.H>

using namespace std;



/** Constructor */

ppFixedU::ppFixedU( const std::string& eqName,
		    const std::string& bdName,
		    const latticeMesh& mesh,	  
		    const scalarField& rho,
		    const scalarField& T,
		    const vectorField& U,
		    pdfField& pdf )
    
    : ppBndCond(mesh, rho, T, U, pdf, bdName, "fixedU") {

    
    // Load boundary value

    dictionary dict("start/boundaries");

    vector<scalar> val = dict.lookUp< vector<scalar> >( eqName + "/" + bdName + "/value" );


    // Resize values at boundary and assign

    _bndVal.resize( _nodes.size() );

    for( uint i = 0 ; i < _bndVal.size() ; i++ ) {

    	_bndVal[i].resize(3);

    	for(uint j = 0 ; j < 3 ; j++)
    	    _bndVal[i][j] = val[j];

    }





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

ppFixedU::~ppFixedU() {}








/** Update pdf field */

void ppFixedU::update( const pseudoPotEquation* ppeq ) {

    
    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();

    vector<scalar> f_eq_nb(q);

    vector<scalar> f_eq_bnd(q);

    


    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {
	

    	uint id = _nodes[i];

    	vector<scalar> Uw = { _bndVal[i][0], _bndVal[i][1], _bndVal[i][2] };

	vector<scalar> nvel = { 0,0,0 };


	// Density and velocity at neighbour node

	const uint nbid = _nbid[i];
	
	scalar nbrho = ppeq->localDensity(nbid);	

	ppeq->localVelocity(nvel, nbid, true);



	// Equilibrium

	ppeq->eqPS( f_eq_nb, nbrho, nvel );

	// ppeq->eqPS( f_eq_bnd, _rho.at(id), Uw );
	ppeq->eqPS( f_eq_bnd, nbrho, Uw );


	
	// Update distribution

    	for( uint k = 0 ; k < q ; k++ ) {	    	    		       		   		    
			
	    _pdf.set(id, k, f_eq_bnd[k] + (_pdf[nbid][k] - f_eq_nb[k] ) );

    	}
	

    }
    
}






// /** Update pdf field */

// void ppFixedU::update( const pseudoPotEquation* ppeq ) {

    
//     // Lattice constants
    
//     const uint q = _mesh.lmodel()->q();

//     vector<scalar> f_eq_nb(q);

//     vector<scalar> f_eq_bnd(q);

//     const vector< vector<int> >& nb = _mesh.nbArray();

//     const vector<uint> reverse = _mesh.lmodel()->reverse();

//     const vector< vector<int> > vel = _mesh.lmodel()->lvel();
    


//     // Move over boundary elements

//     for( uint i = 0 ; i < _nodes.size() ; i++ ) {
	

//     	uint id = _nodes[i];

//     	vector<scalar> Uw = { _bndVal[i][0], _bndVal[i][1], _bndVal[i][2] };

// 	vector<scalar> nvel = { 0,0,0 };


//     	for( uint k = 0 ; k < q ; k++ ) {


//     	    if ( nb[id][k] == -1 ) {

		
//     		// Need density and velocity at neighbour (reverse) node
		    
//     		int nbid = nb[id][ reverse[k] ];


//     		if( nbid != -1 ) {


//     		    // Density and velocity at neighbour node
		    
//     		    scalar nbrho = ppeq->localDensity(nbid);

// 		    ppeq->localVelocity(nvel, nbid, true);

	    
	    			
//     		    // Equilibrium

//     		    ppeq->eqPS( f_eq_nb, nbrho, nvel );

//     		    ppeq->eqPS( f_eq_bnd, nbrho, Uw );
		    
		    

//     		    // Update distribution
			
//     		    _pdf.set(id, k, f_eq_bnd[k] + (_pdf[nbid][k] - f_eq_nb[k] ) );

//     		}
		

//     	    }



//             // Correction for velocities lying on face
	    
// 	    else {

// 		if( nb[id][reverse[k]] != -1 ) {
		
// 		    if( q == 9 ) {

// 			if( ( k == 1 )  ||  ( k == 2 )  ||  ( k == 3 )  ||  ( k == 4 )  ) {

// 			    _pdf[id][k] = 0.5 * _pdf[id][k]  +  0.5 * _pdf[id][reverse[k]];
		
// 			    _pdf[id][reverse[k]] = _pdf[id][k];

// 			}

// 		    }



// 		    else {

// 		    	if( q == 15 ) {

// 		    	    if( ( k == 1 )  ||  ( k == 2 )  ||  ( k == 3 )  ||  ( k == 4 )   ||  ( k == 5 )   ||  ( k == 6 ) ) {

// 		    		_pdf[id][k] = 0.5 * _pdf[id][k]  +  0.5 * _pdf[id][reverse[k]];
		
// 		    		_pdf[id][reverse[k]] = _pdf[id][k];

// 		    	    }

// 		    	}

// 		    }

		    

// 		}

// 	    }
	    
	    

//     	}

//     }
    
// }
