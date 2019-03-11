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





    dictionary dict("properties/macroProperties");

    _Gads = dict.lookUpOrDefault<scalar>( eqName + "/Forces/Adhesive/Gads", 0 );    
    


    // // Read constants for each boundary

    // const map< string, vector<uint> >& boundary = _mesh.boundaries();

    // dictionary dict("start/boundaries");

    // map<string, scalar> Gads;
        
    // for( const auto &bd : boundary ) {
	
    // 	scalar g_ads = dict.lookUpOrDefault<scalar>( eqName + "/" + bd.first + "/Gads", 0 );

    // 	if(  ( g_ads != 0 )  &&  ( bd.second.size() != 0 )  )
    // 	    Gads[bd.first] = g_ads;

    // }



    // // Create constants for closestNodes map

    // const map<string, vector<uint>>& closest = _mesh.nodesCloseToBoundary();
    
    // for(const auto &g : Gads) {

    // 	if( closest.find(g.first) != closest.end() ) {

    // 	    for( auto id : closest.at(g.first) ) {

    // 		closestNodes[id] = g.second;

    // 	    }

    // 	}

    // }

	

}    


/** Destructor */

simpleAdhesive::~simpleAdhesive() {}



/** Interaction potential */

const scalar simpleAdhesive::potential( const scalar& rho, const scalar& T, const scalar& cs2 ) const {

    scalar a = 2 * (eos->p_eos(rho,T) - rho * cs2)  /  (-1);

    scalar b(0);

    (a >= 0)  ?	 b = sqrt(a)  :  b = sqrt(-a);
    
    return b;    

}



/** Force at specific node */

const void simpleAdhesive::force( const uint& i, const scalarField& rho, const scalarField& T,  scalar f[3] ) const {


    f[0] = 0;
    f[1] = 0;
    f[2] = 0;
       

	
    // Reference to neighbour array

    const vector< vector<int> >& nb = _mesh.nbArray();

    

    // Lattice model properties

    const vector< vector<int> >& vel = _mesh.lmodel()->lvel();

    const vector<uint>& reverse = _mesh.lmodel()->reverse();

    const scalar cs2 = _mesh.lmodel()->cs2();

    const uint q = _mesh.lmodel()->q();


    scalar _Fads[3] = {0, 0, 0};
 
    scalar phi = potential( rho.at(i), T.at(i), cs2 );
    

    for( uint k = 1 ; k < q ; k++ ) {
       
    	int neighId = nb[i][ reverse[k] ];

    	if(neighId == -1) {
	    
    	    for( uint j = 0 ; j < 3 ; j++ ) {

    		_Fads[j] +=  _weights[k] * vel[k][j] ;

    	    }

    	}
    

    }

    

    // Extra constant		
    
    for( uint j = 0 ; j < 3 ; j++)
    	f[j] =  -_Gads * phi * phi * _Fads[j] ;   
	

}
