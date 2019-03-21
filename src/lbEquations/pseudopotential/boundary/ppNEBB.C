#include <ppNEBB.H>

using namespace std;



/** Constructor */

ppNEBB::ppNEBB( const std::string& eqName,
		const std::string& bdName,
		const latticeMesh& mesh,	  
		const scalarField& rho,
		const scalarField& T,
		const vectorField& U,
		pdfField& pdf )
    
    : ppWetNodeBnd(eqName, bdName, "NEBB", mesh, rho, T, U, pdf) {

    
}


/** Destructor */

ppNEBB::~ppNEBB() {}




/** Update pdf field */

void ppNEBB::update( const pseudoPotEquation* ppeq ) {


    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();  

    const vector< vector<int> >& nb = _mesh.nbArray();

    const vector<uint> reverse = _mesh.lmodel()->reverse();



    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {


	vector<scalar> Uw = { _bndVal[i][0], _bndVal[i][1], _bndVal[i][2] };

    	uint id = _nodes[i];
	
	scalar Ft[3] = {0,0,0};
	

	// Hand coded

	switch( _mesh.lmodel()->type() ) {

	case latticeModel::latticeType::D2Q9:


	    // Total force at node
	
	    ppeq->totalForce(Ft, id);
	

	    switch( _normal[i] ) {

	    case latticeMesh::normalType::Y1:
	    
		_pdf[id][4] = _pdf[id][2];

		_pdf[id][7] = _pdf[id][5] + 0.5 * (  _pdf[id][1] - _pdf[id][3] ) + 0.25*(Ft[0]+Ft[1]);

		_pdf[id][8] = _pdf[id][6] - 0.5 * (  _pdf[id][1] - _pdf[id][3] ) - 0.25*(Ft[0]-Ft[1]);		
	    
		break;


	    case latticeMesh::normalType::Y0:
	    	
		_pdf[id][2] = _pdf[id][4];

		_pdf[id][5] = _pdf[id][7] - 0.5 * (  _pdf[id][1] - _pdf[id][3] ) - 0.25*(Ft[0]+Ft[1]);

		_pdf[id][6] = _pdf[id][8] + 0.5 * (  _pdf[id][1] - _pdf[id][3] ) + 0.25*(Ft[0]-Ft[1]);
				
		break;


	    case latticeMesh::normalType::X0:

		_pdf[id][1] = _pdf[id][3];

		_pdf[id][5] = _pdf[id][7] + 0.5 * (  _pdf[id][4] - _pdf[id][2] ) + 0.25*(Ft[0]+Ft[1]);

		_pdf[id][8] = _pdf[id][6] - 0.5 * (  _pdf[id][4] - _pdf[id][2] ) + 0.25*(Ft[0]-Ft[1]);
			    
	    	break;


	    case latticeMesh::normalType::X1:

		_pdf[id][3] = _pdf[id][1];

		_pdf[id][7] = _pdf[id][5] - 0.5 * (  _pdf[id][4] - _pdf[id][2] ) - 0.25*(Ft[0]+Ft[1]);

		_pdf[id][6] = _pdf[id][8] + 0.5 * (  _pdf[id][4] - _pdf[id][2] ) - 0.25*(Ft[0]-Ft[1]);
		    
	    	break;	    
	    

	    default:

		// Apply simple bounce back rule

		for(uint k = 0 ; k < q ; k++) {

		    if ( nb[id][k] == -1 ) {

			_pdf[id][k] = _pdf[id][reverse[k]];

		    }		    

		}
	    
		break;

	    }

	
	    break;


	    

	default:

	    cout << "Non-equilibrium bounce-back not yet implemented for this lattice type" << endl;

	    exit(1);

	    break;
	
	}
	

    }
    
}
