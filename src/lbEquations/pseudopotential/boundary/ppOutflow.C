#include <ppOutflow.H>

using namespace std;


/** Constructor */

ppOutflow::ppOutflow( const string& eqName,
		      const string& bdName,
		      const latticeMesh& mesh,
		      const scalarField& rho,
		      const scalarField& T,
		      const vectorField& U,
		      pdfField& pdf )

    : ppFixedU( eqName, bdName, mesh, rho, T, U, pdf ) {



    // Resize neighbour indices, compute normals and check indexing

    _nbid.resize( _bndVal.size() );

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

ppOutflow::~ppOutflow() {}


/** Update pdf field */

void ppOutflow::update( const pseudoPotEquation* ppeq ) {

    
    // First compute value over boundary according to _grad

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {

	for( uint j = 0 ; j < 3 ; j++ ) {
	
	    _bndVal[i][j] = _U.at(_nodes[i],j);

	}

    }


    // Compute base version
    
    ppFixedU::update(ppeq);


}
