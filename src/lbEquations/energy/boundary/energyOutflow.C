#include <energyOutflow.H>

using namespace std;


/** Constructor */

energyOutflow::energyOutflow( const string& eqName,
		      const string& bdName,
		      const latticeMesh& mesh,
		      const scalarField& rho,
		      const scalarField& T,
		      const vectorField& U,
		      pdfField& pdf )

    : energyBndCond(mesh, rho, T, U, pdf, bdName) {


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

energyOutflow::~energyOutflow() {}




/** Update pdf field */

void energyOutflow::update( const energyEquation* eeq ) {


    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();

    // const vector< vector<int> >& nb = _mesh.nbArray();

    vector<scalar> nvel = { 0,0,0 };


    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {

	
	uint id = _nodes[i];

	uint nid = _nbid[i];
		
	
	// Velocity at neighbour node		    

	scalar uAdv(0);
	
	for(uint j = 0 ; j < 3 ; j++)
	    uAdv += _normal[i][j] * _U.at(nid)[j];


	// Update unknowun distributions for f
	    
	for( uint k = 0 ; k < q ; k++ ) {
		
	    _pdf.set( id,
	    	      k,
	    	      ( _pdf[id][k] + uAdv*_pdf[nid][k] ) / (1+uAdv)
	    	);	    

	}

	
    }

    

}
