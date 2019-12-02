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






    // Load mean normal

    dictionary dict("start/boundaries");

    _normal = dict.lookUp< vector<scalar> >( eqName + "/" + bdName + "/meanNormal" );

    

    // Detect correspondence with lattice velocities

    ln = 0;

    const uint q = _mesh.lmodel()->q();

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();    

    for( uint k = 0 ; k < q ; k++ ) {
	    
	if(          ( vel[k][0] == (int)_normal[0] )
		 &&  ( vel[k][1] == (int)_normal[1] )
		 &&  ( vel[k][2] == (int)_normal[2] )  ) {
		
	    ln = k;

	}

    }

    if(ln == 0) {
	
	cout << " [ERROR]  Unable to detect normal on nodes over boundary " << bdName;

	exit(1);

    }



    // Pre-allocate density map

    for( auto n : _nodes )
	_rhow[n] = _rho.at(n);
    

    

}


/** Destructor */

ppOutflowWithNEBB::~ppOutflowWithNEBB() {}



/** Update pdf field */

void ppOutflowWithNEBB::update( const pseudoPotEquation* ppeq ) {


    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();

    const vector< vector<int> >& nb = _mesh.nbArray();

    vector<scalar> nvel = { 0,0,0 };


    
    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {

	
    	uint id = _nodes[i];

    	uint nid = nb[id][ln];
		
	
    	// Velocity at neighbour node		    

    	ppeq->localVelocity(nvel, nid, true);

    	scalar uAdv(0);
	
    	for(uint j = 0 ; j < 3 ; j++)
    	    uAdv += _normal[j] * nvel[j];



	// Transport density to boundary

	_rhow[id] = ( _rho.at(id) + uAdv*ppeq->localDensity(nid) ) / (1+uAdv);


	// Transport velocity to boundary

	for(uint j = 0 ; j < 3 ; j++)
	    _bndVal[i][j] = ( _U.at(id)[j] + uAdv*nvel[j] ) / (1+uAdv);

	
    }



    
    // Apply general NEBB

    ppGeneralNEBB::update( ppeq, _rhow );

}
