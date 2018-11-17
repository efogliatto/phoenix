#include <simpleAdhesive.H>

#include <dictionary.H>


using namespace std;


/** Constructor */

simpleAdhesive::simpleAdhesive( const string& dictName,
				const string& eqName,
				const latticeMesh& mesh ) 

    
    : adhesiveForce(dictName, eqName, mesh) {



    // Update pseudo pot weights

    string model = mesh.lmodel()->name();

    if( model == "D2Q9" ) {

    	_weights.push_back( 0 );
    	_weights.push_back( 1.0/9 );
    	_weights.push_back( 1.0/9 );
    	_weights.push_back( 1.0/9 );
    	_weights.push_back( 1.0/9 );
    	_weights.push_back( 1.0/36 );
    	_weights.push_back( 1.0/36 );
    	_weights.push_back( 1.0/36 );
    	_weights.push_back( 1.0/36 );

    }

    else {

	if( model == "D3Q15" ) {

	    _weights.push_back( 0 );

	    for( uint i = 1 ; i <= 6 ; i++ )
		_weights.push_back( 1.0 / 9.0 );

	    for( uint i = 7 ; i < 15 ; i++ )
		_weights.push_back( 1.0 / 72.0 );
	    


	}

	else {

	    cout << " [ERROR]  Potential weights not implemented yet for " << model << endl;

	}

    }





    // Read constants for each boundary

    const map< string, vector<uint> >& boundary = _mesh.boundaries();

    dictionary dict("start/boundaries");

    
    // // Lattice model properties

    // const vector< vector<int> >& nb = _mesh.nbArray();
	    
    // const uint q = _mesh.lmodel()->q();    


    
    for( const auto &bd : boundary ) {
	
	scalar g_ads = dict.lookUpOrDefault<scalar>( eqName + "/" + bd.first + "/Gads", 0 );

	
	// // Check for closest nodes. Move over boundary nodes and check if neighbours are on boundaries
	
	// if(  ( g_ads != 0 )  ) {	   
	  
	    
	//     for( const auto id : bd.second ) {		
		
	// 	for(uint k = 1 ; k < q ; k++) {
	    
	// 	    int aux = nb[id][k];

	// 	    if(aux != -1) {

	// 		bool is_on_bnd(false);

	// 		for(uint l = 1 ; l < q ; l++) {

	// 	    	    if( nb[aux][l] == -1 ) {
				
	// 	    		is_on_bnd = true;

	// 	    	    }

	// 		}

	// 		if( (!is_on_bnd)  &&  (closestNodes.find(id) == closestNodes.end()) ) {

	// 		    closestNodes[id] = g_ads;
			    
	// 		}
			

	// 	    }
		    
		    
		    
		    

	// 	}
		

	//     }
	    

	// }

    }

}    


/** Destructor */

simpleAdhesive::~simpleAdhesive() {}



/** Interaction potential */

const scalar simpleAdhesive::potential( const scalar& rho, const scalar& T, const scalar& cs2 ) const {

    scalar a = 2 * (eos->p_eos(rho,T) - rho * cs2);

    scalar b(0);

    (a >= 0)  ?	 b = sqrt(a)  :  b = sqrt(-a);
    
    return b;

}



/** Force at specific node */

const void simpleAdhesive::force( const uint& i, const scalarField& rho, const scalarField& T,  scalar f[3] ) const {

    
    // if( closestNodes.find(i) != closestNodes.end() ) {

	
    // 	// Reference to neighbour array

    // 	const vector< vector<int> >& nb = _mesh.nbArray();


    // 	// Lattice model properties

    // 	const vector< vector<int> >& vel = _mesh.lmodel()->lvel();

    // 	const vector<uint>& reverse = _mesh.lmodel()->reverse();

    // 	const scalar cs2 = _mesh.lmodel()->cs2();

    // 	const uint q = _mesh.lmodel()->q();


    // 	vector<scalar> F = {0, 0, 0};	    

    // 	for( uint k = 1 ; k < q ; k++ ) {
       
    // 	    int neighId = nb[i][ reverse[k] ];

    // 	    scalar _rho = rho.at(neighId);

    // 	    scalar _T = T.at(neighId);
	    
    // 	    scalar alpha = _weights[k] * potential( _rho, _T, cs2 );

	    
    // 	    for( uint j = 0 ; j < 3 ; j++ ) {

    // 		F[j] +=  alpha * (scalar)vel[k][j] ;

    // 	    }
    

    // 	}

		

	// // Extra constant
		
	// scalar beta = potential( rho.at(i), T.at(i), cs2 ) * closestNodes.at(i);
	    
    
	// for( uint j = 0 ; j < 3 ; j++) {
	
	//     f[j] =  F[j] * beta;
	
	// }	
	

    // }

}
