#include <ppGeneralNEBB.H>

#include <math.h>

using namespace std;



// Rotation matrix along unitary vector

const vector< vector<scalar> > ppGeneralNEBB::rotationMatrix( const scalar u, const scalar v, const scalar w, const scalar theta ) const {

    vector< vector<scalar> > _rot = {{0,0,0},{0,0,0},{0,0,0}};

    scalar ct( cos(theta) );
    scalar st( sin(theta) );
    
    _rot[0][0] = ct + u*u*(1-ct);
    _rot[0][1] = u*v*(1-ct) - w*st;
    _rot[0][2] = u*w*(1-ct) + v*st;
    _rot[1][0] = u*v*(1-ct) + w*st;
    _rot[1][1] = ct + v*v*(1-ct);
    _rot[1][2] = v*w*(1-ct)-u*st;
    _rot[2][0] = u*w*(1-ct) - v*st;
    _rot[2][1] = v*w*(1-ct) + u*st;
    _rot[2][2] = ct + w*w*(1-ct);
    
    return _rot;

}



/** Constructor */

ppGeneralNEBB::ppGeneralNEBB( const std::string& eqName,
			      const std::string& bdName,
			      const latticeMesh& mesh,	  
			      const scalarField& rho,
			      const scalarField& T,
			      const vectorField& U,
			      pdfField& pdf )
    
    : ppWetNodeBnd(eqName, bdName, "generalNEBB", mesh, rho, T, U, pdf) {


    // Allocate lattice transformation indices

    switch( _mesh.lmodel()->type() ) {

    case latticeModel::latticeType::D2Q9:

	rotation[ latticeMesh::normalType::X0 ] = rotationMatrix(0,0,1,M_PI/2);
	rotation[ latticeMesh::normalType::X1 ] = rotationMatrix(0,0,1,-M_PI/2);
	rotation[ latticeMesh::normalType::Y0 ] = rotationMatrix(0,0,0,0);
	rotation[ latticeMesh::normalType::Y1 ] = rotationMatrix(0,0,1,M_PI);
	
	// trIndex[ latticeMesh::normalType::Y0 ] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
	// trIndex[ latticeMesh::normalType::Y1 ] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
	// trIndex[ latticeMesh::normalType::X0 ] = { 0, 4, 1, 2, 3, 8, 5, 6, 7 };
	// trIndex[ latticeMesh::normalType::X1 ] = { 0, 2, 3, 4, 1, 6, 7, 8, 5 };
	
	break;


    case latticeModel::latticeType::D3Q15:

	rotation[ latticeMesh::normalType::X0 ] = rotationMatrix(0, 0, 0, 0);
	rotation[ latticeMesh::normalType::X1 ] = rotationMatrix(0, 1, 0, M_PI);
	rotation[ latticeMesh::normalType::Y0 ] = rotationMatrix(0, 0, 1,-M_PI/2);
	rotation[ latticeMesh::normalType::Y1 ] = rotationMatrix(0, 0, 1, M_PI/2);
	rotation[ latticeMesh::normalType::Z0 ] = rotationMatrix(0, 1, 0, M_PI/2);
	rotation[ latticeMesh::normalType::Z1 ] = rotationMatrix(0, 1, 0,-M_PI/2);	

	break;

	
    default:

	cout << " [ERROR] Unable to use lattice transformation coefficients in generalNEBB" << endl << endl;

	exit(1);

	break;

    }




    // Move over lattice velocities, transform and check velocity indexing

    const vector< vector<int> >& vel = _mesh.lmodel()->lvel();


    for( auto const& rot : rotation ) {
	

    	trIndex[rot.first].resize( _mesh.lmodel()->q(),0 );

	
    	for( uint k = 0 ; k < _mesh.lmodel()->q() ; k++ ) {


    	    // Rotate grid velocity
	    
    	    scalar rotVel[3] = {0,0,0};

    	    for( uint i = 0 ; i < 3 ; i++ ) {

    		for( uint j = 0 ; j < 3 ; j++ ) {

    		    rotVel[i] = rotVel[i] + rot.second[i][j]*vel[k][j];

    		}


		// Round for better comparison
		
		rotVel[i] = std::round( rotVel[i] );

		if( rotVel[i] == -0 )
		    rotVel[i] = 0;

    	    }
	   	   



    	    // Check for grid correspondence

    	    for( uint i = 0 ; i < _mesh.lmodel()->q() ; i++ ) {

    		if( vel[i][0] == (int)rotVel[0] ) {
		    
		    if( vel[i][1] == (int)rotVel[1] ) {

			if ( vel[i][2] == (int)rotVel[2] ) {
			
			    trIndex[ rot.first ][i] = k;

			}

		    }
		    
    		}

    	    }


	    

    	}

    }

    
    
    
}


/** Destructor */

ppGeneralNEBB::~ppGeneralNEBB() {}




/** Update pdf field */

void ppGeneralNEBB::update( const pseudoPotEquation* ppeq ) {


    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();  

    const vector< vector<int> >& nb = _mesh.nbArray();

    const vector<uint> reverse = _mesh.lmodel()->reverse();



    // Auxiliary arrays

    vector<scalar> Uw = {0,0,0};	
	
    scalar Ft[3] = {0,0,0};

    scalar rotFt[3] = {0,0,0};

        

    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {


	// Boundary node index (on lattice indexing)
	
	uint id = _nodes[i];

	

	// Apply simple bounce back rule for corners
	
	if( _normal[i] == latticeMesh::normalType::UNDEF ) {

	    for(uint k = 0 ; k < q ; k++) {

		if ( nb[id][k] == -1 ) {

		    _pdf[id][k] = _pdf[id][reverse[k]];

		}		    

	    }

	}



	else {
	    
	    // Total force at node
	
	    ppeq->totalForce(Ft, id);

	    // for(uint j = 0 ; j < 3 ; j++)	    
	    // 	Uw[j] = _bndVal[i][j];


	    // Rotate force and velocity at node

	    vector< vector<scalar> >& rot = rotation[ _normal[i] ];	    

	    for(uint j = 0 ; j < 3 ; j++)	    
	    	Uw[j] = rot[j][0] * _bndVal[i][0]  +  rot[j][1] * _bndVal[i][1]  +  rot[j][2] * _bndVal[i][2];

	    for(uint j = 0 ; j < 3 ; j++)	    
	    	rotFt[j] = rot[j][0] * Ft[0]  +  rot[j][1] * Ft[1]  +  rot[j][2] * Ft[2];	    

	    

	    
	    

	    // Use transformation vectors for general formula
	    
	    const vector<uint>& tr = trIndex[ _normal[i] ];

	    scalar rho_w(0);

	    scalar Delta[3] = {0,0,0};



	    switch( _mesh.lmodel()->type() ) {
		
	    
	    case latticeModel::latticeType::D2Q9:


		// Wall density
		
		rho_w = ( _pdf[id][tr[0]] + _pdf[id][tr[1]] + _pdf[id][tr[3]] + 2.0*(_pdf[id][tr[4]] + _pdf[id][tr[7]] + _pdf[id][tr[8]]) - 0.5*rotFt[1]  ) / (1.0 - Uw[1]);


		// Update distribution
		
		_pdf[id][ tr[2] ] = _pdf[id][ tr[4] ]  +  (2.0/3.0) * rho_w  *  Uw[1];

		_pdf[id][ tr[5] ] = _pdf[id][ tr[7] ]
		                  + 0.5 * rho_w * Uw[0]
		                  + (rho_w/6.0) * Uw[1]
		                  - 0.5 * ( _pdf[id][ tr[1] ] - _pdf[id][ tr[3] ] )
		                  - 0.25 * ( rotFt[0] + rotFt[1] );

		_pdf[id][ tr[6] ] = _pdf[id][ tr[8] ]
		                  - 0.5 * rho_w * Uw[0]
		                  + (rho_w/6.0) * Uw[1]
		                  + 0.5 * ( _pdf[id][ tr[1] ] - _pdf[id][ tr[3] ] )
		                  - 0.25 * ( -rotFt[0] + rotFt[1] );		


		break;





	    case latticeModel::latticeType::D3Q15:


		// Wall density
		
		rho_w = ( _pdf[id][tr[0]] + _pdf[id][tr[3]] + _pdf[id][tr[4]] + _pdf[id][tr[5]] + _pdf[id][tr[6]]  +  2*( _pdf[id][tr[2]] + _pdf[id][tr[8]] + _pdf[id][tr[10]] + _pdf[id][tr[12]] + _pdf[id][tr[14]] )  - 0.5*rotFt[0] )  /  (1 - Uw[0]);


		// Momentum correctors
		
		Delta[0] = -rotFt[0] / 8.0;

		Delta[1] = -(_pdf[id][tr[3]] - _pdf[id][tr[4]]) / 4.0   +   rho_w * Uw[1] / 6.0   -   rotFt[1] / 8.0;

		Delta[2] = -(_pdf[id][tr[5]] - _pdf[id][tr[6]]) / 4.0   +   rho_w * Uw[2] / 6.0   -   rotFt[2] / 8.0;


		// Update distribution

		_pdf[id][tr[1]]  = _pdf[id][tr[2]]     +    2.0 * rho_w * Uw[0] / 3.0;

		_pdf[id][tr[7]]  = _pdf[id][tr[14]]    +    rho_w * ( Uw[0] + Uw[1] + Uw[2] ) / 12.0    +    Delta[0] + Delta[1] + Delta[2];

		_pdf[id][tr[9]]  = _pdf[id][tr[12]]    +    rho_w * ( Uw[0] - Uw[1] + Uw[2] ) / 12.0    +    Delta[0] - Delta[1] + Delta[2];

		_pdf[id][tr[11]] = _pdf[id][tr[10]]    +    rho_w * ( Uw[0] + Uw[1] - Uw[2] ) / 12.0    +    Delta[0] + Delta[1] - Delta[2];

		_pdf[id][tr[13]] = _pdf[id][tr[8]]     +    rho_w * ( Uw[0] - Uw[1] - Uw[2] ) / 12.0    +    Delta[0] - Delta[1] - Delta[2];
		
		

		break;


	    }

	}
	

    }

    
}
