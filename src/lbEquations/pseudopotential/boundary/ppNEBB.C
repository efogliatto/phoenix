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

    vector<scalar> f_eq(q);

    


    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {


	vector<scalar> Uw = { _bndVal[i][0], _bndVal[i][1], _bndVal[i][2] };

    	uint id = _nodes[i];

	scalar rhow(0);
	
	scalar Ft[3] = {0,0,0};
	

	// Hand coded

	switch( _mesh.lmodel()->type() ) {

	case latticeModel::latticeType::D2Q9:


	    // Total force at node
	
	    ppeq->totalForce(Ft, id);
	

	    switch( _normal[i] ) {

	    case latticeMesh::normalType::Y1:
	    
		rhow = ppeq->localDensityWithUnknowns( id, _normal[i] );

		ppeq->eqPS( f_eq, rhow, Uw );

		_pdf[id][4] = _pdf[id][2] + (f_eq[4] - f_eq[2]);

		_pdf[id][7] = _pdf[id][5] + 0.5 * (  -rhow*Uw[0] +  _pdf[id][1] - _pdf[id][3] + f_eq[7] - f_eq[5] + f_eq[8] - f_eq[6] + 0.5*Ft[0]  );

		_pdf[id][8] = _pdf[id][6] + 0.5 * (   rhow*Uw[0] -  _pdf[id][1] + _pdf[id][3] + f_eq[7] - f_eq[5] + f_eq[8] - f_eq[6] - 0.5*Ft[0]  );
	    
		break;


	    case latticeMesh::normalType::Y0:
	    
		rhow = ppeq->localDensityWithUnknowns( id, _normal[i] );

		ppeq->eqPS( f_eq, rhow, Uw );

		_pdf[id][2] = _pdf[id][4] + (f_eq[2] - f_eq[4]);

		_pdf[id][5] = _pdf[id][7] + 0.5 * (   rhow*Uw[0] -  _pdf[id][1] + _pdf[id][3] - f_eq[7] + f_eq[5] - f_eq[8] + f_eq[6] - 0.5*Ft[0]  );

		_pdf[id][6] = _pdf[id][8] + 0.5 * (  -rhow*Uw[0] +  _pdf[id][1] - _pdf[id][3] - f_eq[7] + f_eq[5] - f_eq[8] + f_eq[6] + 0.5*Ft[0]  );
		
		break;


	    case latticeMesh::normalType::X0:

		rhow = ppeq->localDensityWithUnknowns( id, _normal[i] );

		ppeq->eqPS( f_eq, rhow, Uw );

		_pdf[id][1] = _pdf[id][3] + (f_eq[1] - f_eq[3]);

		_pdf[id][5] = _pdf[id][7] + 0.5 * (   rhow*Uw[1] -  _pdf[id][2] + _pdf[id][4] + f_eq[5] - f_eq[7] + f_eq[8] - f_eq[6] - 0.5*Ft[1]  );

		_pdf[id][8] = _pdf[id][6] + 0.5 * (  -rhow*Uw[1] +  _pdf[id][2] - _pdf[id][4] + f_eq[5] - f_eq[7] + f_eq[8] - f_eq[6] + 0.5*Ft[1]  );		
	    
	    	break;


	    case latticeMesh::normalType::X1:

		rhow = ppeq->localDensityWithUnknowns( id, _normal[i] );

		ppeq->eqPS( f_eq, rhow, Uw );

		_pdf[id][3] = _pdf[id][1] + (f_eq[3] - f_eq[1]);

		_pdf[id][7] = _pdf[id][5] + 0.5 * (  -rhow*Uw[1] +  _pdf[id][2] - _pdf[id][4] - f_eq[5] + f_eq[7] - f_eq[8] + f_eq[6] + 0.5*Ft[1]  );

		_pdf[id][6] = _pdf[id][8] + 0.5 * (   rhow*Uw[1] -  _pdf[id][2] + _pdf[id][4] - f_eq[5] + f_eq[7] - f_eq[8] + f_eq[6] - 0.5*Ft[1]  );
	    
	    	break;	    
	    

	    default:

		cout << "Unable to compute local density" << endl;

		exit(1);
	    
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
