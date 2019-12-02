#include <ppOutflowWithNEBB.H>

using namespace std;


/** Constructor */

ppOutflowWithNEBB::ppOutflowWithNEBB( const string& eqName,
				      const string& bdName,
				      const latticeMesh& mesh,
				      const scalarField& rho,
				      const scalarField& T,
				      const vectorField& U,
				      pdfField& pdf )

    : ppGeneralNEBB(eqName, bdName, mesh, rho, T, U, pdf) {


    // Resize neighbour indices, compute normals and check indexing

    _nbid.resize( _nodes.size() );

    // const vector< vector<int> >& nb = mesh.nbArray();

    // const vector<uint>& reverse = mesh.lmodel()->reverse();

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();

    const uint q = mesh.lmodel()->q();



    // Load mean normal

    dictionary dict("start/boundaries");

    vector<scalar> normal = dict.lookUp< vector<scalar> >( eqName + "/" + bdName + "/meanNormal" );


    // Detect correspondence with lattice velocities

    ln = 0;

    for( uint k = 0 ; k < q ; k++ ) {
	    
	if(          ( vel[k][0] == (int)normal[0] )
		 &&  ( vel[k][1] == (int)normal[1] )
		 &&  ( vel[k][2] == (int)normal[2] )  ) {
		
	    ln = k;

	}

    }

    if(ln == 0) {
	
	cout << " [ERROR]  Unable to detect normal on nodes over boundary " << bdName;

	exit(1);

    }
    

    

}


/** Destructor */

ppOutflowWithNEBB::~ppOutflowWithNEBB() {}


/** Update pdf field */

void ppOutflowWithNEBB::update( const pseudoPotEquation* ppeq ) {


    // // Lattice constants
    
    // const uint q = _mesh.lmodel()->q();

    // // const vector< vector<int> >& nb = _mesh.nbArray();

    // vector<scalar> nvel = { 0,0,0 };


    // // Move over boundary elements

    // for( uint i = 0 ; i < _nodes.size() ; i++ ) {

	
    // 	uint id = _nodes[i];

    // 	uint nid = _nbid[i];
		
	
    // 	// Velocity at neighbour node		    

    // 	ppeq->localVelocity(nvel, nid, true);

    // 	scalar uAdv(0);
	
    // 	for(uint j = 0 ; j < 3 ; j++)
    // 	    uAdv += _normal[i][j] * nvel[j];



    // 	// Update unknowun distributions for f
	    
    // 	for( uint k = 0 ; k < q ; k++ ) {
		
    // 	    _pdf.set( id,
    // 	  	      k,
    // 		      ( _oldPdf[i][k] + uAdv*_pdf[nid][k] ) / (1+uAdv)
    // 		);

    // 	}

    // 	for( uint k = 0 ; k < q ; k++ )		
    // 	    _oldPdf[i][k] = _pdf[id][k];


	
    // }



    
    // Apply general NEBB

    ppGeneralNEBB::update( ppeq );

}





// /** Update interaction force */
    
// const void ppOutflow::updateIntForce( pseudoPotEquation* ppeq ) const {


//     // Lattice constants    

//     vector<scalar> nvel = { 0,0,0 };

//     vector<scalar> Fint = {0,0,0};


//     // Move over boundary elements

//     for( uint i = 0 ; i < _nodes.size() ; i++ ) {

	
// 	uint id = _nodes[i];

// 	uint nid = _nbid[i];
		
	
// 	// Velocity at neighbour node		    

// 	ppeq->localVelocity(nvel, nid, true);

// 	scalar uAdv(0);
	
// 	for(uint j = 0 ; j < 3 ; j++)
// 	    uAdv += _normal[i][j] * nvel[j];



// 	// Update 
	    
// 	for( uint j = 0 ; j < 3 ; j++ )
// 	    Fint[j] = ppeq->intForce(id,j)   +   uAdv * ppeq->intForce(nid,j) / (1+uAdv);

	
// 	ppeq->setIntForce( id, Fint );
	
	    
	
//     }

    

// }
