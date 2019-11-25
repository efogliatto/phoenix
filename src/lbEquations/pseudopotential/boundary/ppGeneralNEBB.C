#include <ppGeneralNEBB.H>

using namespace std;



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

	trIndex[ latticeMesh::normalType::Y0 ] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
	trIndex[ latticeMesh::normalType::Y1 ] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
	trIndex[ latticeMesh::normalType::X0 ] = { 0, 4, 1, 2, 3, 8, 5, 6, 7 };
	trIndex[ latticeMesh::normalType::X1 ] = { 0, 2, 3, 4, 1, 6, 7, 8, 5 };

	prDirIndex[ latticeMesh::normalType::Y0 ] = {  1, 2 };
	prDirIndex[ latticeMesh::normalType::Y1 ] = {  1, 2 };
	prDirIndex[ latticeMesh::normalType::X0 ] = {  2, 1 };
	prDirIndex[ latticeMesh::normalType::X1 ] = {  2, 1 };

	prDirSign[ latticeMesh::normalType::Y0 ] = {  1,  1 };
	prDirSign[ latticeMesh::normalType::Y1 ] = { -1, -1 };
	prDirSign[ latticeMesh::normalType::X0 ] = { -1,  1 };
	prDirSign[ latticeMesh::normalType::X1 ] = {  1, -1 };	
	
	break;


    case latticeModel::latticeType::D3Q15:



	break;

	
    default:

	cout << " [ERROR] Unable to use lattice transformation coefficients in generalNEBB" << endl << endl;

	exit(1);

	break;

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
	    
	    // Total force  and velocity at node
	
	    ppeq->totalForce(Ft, id);

	    for(uint j = 0 ; j < 3 ; j++)	    
		Uw[j] = _bndVal[i][j];

	    
	    

	    // Use transformation vectors for general formula
	    
	    const vector<uint>& tr = trIndex[ _normal[i] ];

	    const vector<int>& pr = prDirIndex[ _normal[i] ];

	    const vector<int>& pr_sign = prDirSign[ _normal[i] ];


	    switch( _mesh.lmodel()->type() ) {
	    
	    case latticeModel::latticeType::D2Q9:


		scalar rho_w = ( _pdf[id][tr[0]] + _pdf[id][tr[1]] + _pdf[id][tr[3]] + 2.0*(_pdf[id][tr[4]] + _pdf[id][tr[7]] + _pdf[id][tr[8]]) - 0.5*pr_sign[1]*Ft[pr[1]]  ) / (1.0 - pr_sign[1]*Uw[pr[1]]);
	    
		_pdf[id][ tr[2] ] = _pdf[id][ tr[4] ]  +  (2.0/3.0) * rho_w  *  Uw[ pr[1] ] * pr_sign[1];

		_pdf[id][ tr[5] ] = _pdf[id][ tr[7] ]
		    + 0.5 * rho_w * pr_sign[0] * Uw[pr[0]]
		    + (rho_w/6.0) * pr_sign[1] *Uw[pr[1]]
		    - 0.5 * ( _pdf[id][ tr[1] ] - _pdf[id][ tr[3] ] )
		    - 0.25 * ( pr_sign[0]*Ft[pr[0]] + pr_sign[1]*Ft[pr[1]] );

		_pdf[id][ tr[6] ] = _pdf[id][ tr[8] ]
		    - 0.5 * rho_w * pr_sign[0] * Uw[pr[0]]
		    + (rho_w/6.0) * pr_sign[1] *Uw[pr[1]]
		    + 0.5 * ( _pdf[id][ tr[1] ] - _pdf[id][ tr[3] ] )
		    - 0.25 * ( -pr_sign[0]*Ft[pr[0]] + pr_sign[1]*Ft[pr[1]] );


		break;


	    }

	}
	

    }



	
	
	

    // 	// Hand coded

    // 	switch( _mesh.lmodel()->type() ) {

    // 	case latticeModel::latticeType::D2Q9:


    // 	    // Total force at node
	
    // 	    ppeq->totalForce(Ft, id);
	

    // 	    switch( _normal[i] ) {

    // 	    case latticeMesh::normalType::Y1:
	    
    // 		_pdf[id][4] = _pdf[id][2];

    // 		_pdf[id][7] = _pdf[id][5] + 0.5 * (  _pdf[id][1] - _pdf[id][3] ) + 0.25*(Ft[0]+Ft[1]);

    // 		_pdf[id][8] = _pdf[id][6] - 0.5 * (  _pdf[id][1] - _pdf[id][3] ) - 0.25*(Ft[0]-Ft[1]);		
	    
    // 		break;


    // 	    case latticeMesh::normalType::Y0:
	    	
    // 		_pdf[id][2] = _pdf[id][4];

    // 		_pdf[id][5] = _pdf[id][7] - 0.5 * (  _pdf[id][1] - _pdf[id][3] ) - 0.25*(Ft[0]+Ft[1]);

    // 		_pdf[id][6] = _pdf[id][8] + 0.5 * (  _pdf[id][1] - _pdf[id][3] ) + 0.25*(Ft[0]-Ft[1]);
				
    // 		break;


    // 	    case latticeMesh::normalType::X0:

    // 		_pdf[id][1] = _pdf[id][3];

    // 		_pdf[id][5] = _pdf[id][7] + 0.5 * (  _pdf[id][4] - _pdf[id][2] ) + 0.25*(Ft[0]+Ft[1]);

    // 		_pdf[id][8] = _pdf[id][6] - 0.5 * (  _pdf[id][4] - _pdf[id][2] ) + 0.25*(Ft[0]-Ft[1]);
			    
    // 	    	break;


    // 	    case latticeMesh::normalType::X1:

    // 		_pdf[id][3] = _pdf[id][1];

    // 		_pdf[id][7] = _pdf[id][5] - 0.5 * (  _pdf[id][4] - _pdf[id][2] ) - 0.25*(Ft[0]+Ft[1]);

    // 		_pdf[id][6] = _pdf[id][8] + 0.5 * (  _pdf[id][4] - _pdf[id][2] ) - 0.25*(Ft[0]-Ft[1]);
		    
    // 	    	break;	    
	    

    // 	    default:

    // 		// Apply simple bounce back rule

    // 		for(uint k = 0 ; k < q ; k++) {

    // 		    if ( nb[id][k] == -1 ) {

    // 			_pdf[id][k] = _pdf[id][reverse[k]];

    // 		    }		    

    // 		}
	    
    // 		break;

    // 	    }

	
    // 	    break;


	    

    // 	default:

    // 	    cout << "Non-equilibrium bounce-back not yet implemented for this lattice type" << endl;

    // 	    exit(1);

    // 	    break;
	
    // 	}
	

    // }
    
}
