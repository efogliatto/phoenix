#include <ppWetNodeBnd.H>

using namespace std;



/** Constructor */

ppWetNodeBnd::ppWetNodeBnd( const std::string& eqName,
			    const std::string& bdName,
			    const std::string& bdType,
			    const latticeMesh& mesh,	  
			    const scalarField& rho,
			    const scalarField& T,
			    const vectorField& U,
			    pdfField& pdf )
    
    : ppBndCond(mesh, rho, T, U, pdf, bdName, bdType) {

    
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

    _normal.resize( _nodes.size() );

    const vector< vector<int> >& nb = mesh.nbArray();

    const vector<uint>& reverse = mesh.lmodel()->reverse();

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();

    const uint q = mesh.lmodel()->q();


    
    for( uint i = 0 ; i < _nodes.size() ; i++ ) {


	// Node local id
	
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


	    // Assign normal type (outward positive convention)

	    if(  ( normal[0] == -1 )  &&  ( normal[1] == 0 )  &&  ( normal[2] == 0 ) ) {

	    	_normal[i] = normalType::X0;

	    }

	    else {

		if(  ( normal[0] == 1 )  &&  ( normal[1] == 0 )  &&  ( normal[2] == 0 ) ) {

		    _normal[i] = normalType::X1;

		}

		else {

		    if(  ( normal[0] == 0 )  &&  ( normal[1] == -1 )  &&  ( normal[2] == 0 ) ) {

			_normal[i] = normalType::Y0;

		    }

		    else {

			if(  ( normal[0] == 0 )  &&  ( normal[1] == 1 )  &&  ( normal[2] == 0 ) ) {

			    _normal[i] = normalType::Y1;

			}

			else {

			    if(  ( normal[0] == 0 )  &&  ( normal[1] == 0 )  &&  ( normal[2] == -1 ) ) {

				_normal[i] = normalType::Z0;

			    }

			    else {

				if(  ( normal[0] == 0 )  &&  ( normal[1] == 0 )  &&  ( normal[2] == 1 ) ) {

				    _normal[i] = normalType::Z1;

				}

				else {

				    _normal[i] = normalType::UNDEF;

				}

			    }			    				

			}			    

		    }			

		}	
		    

	    }


	    
	    
	    
    	}

    	else {

    	    cout << " [ERROR]  Unable to detect normal on node " << i << " over boundary " << bdName;

    	    exit(1);

    	}
             


	
    	// Detect correspondence with lattice velocities and add normal neighbour

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

ppWetNodeBnd::~ppWetNodeBnd() {}
