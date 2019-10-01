#include <flatCavity.H>

#include <algorithm>


using namespace std;


/** Constructor */
    
flatCavity::flatCavity( const string& bdname ) : cavityModel(bdname) {}


/** Destructor */

flatCavity::~flatCavity() {}


/** Compute adhesive coefficients */

const vector< pair<uint, scalar> > flatCavity::coeffs(const vector<uint>& nodes,
						      const vector< vector<uint> >& loc,
						      const vector< pair<uint,uint> >& spots) const {


    // Coefficients array
    
    vector< pair<uint, scalar> > Gads( nodes.size() );

    for( uint i = 0 ; i < nodes.size() ; i++ )
	Gads[i] = make_pair( nodes[i], max_gads );
    

    
    // Move over spots

    for( auto sp : spots ) {


    	// Move over nodes until spot is found

    	bool find(false);

    	uint i(0);


    	while( (find==false)  &&  (i<nodes.size()) ) {
	    
    	    if( sp.first == nodes[i] ) {
		
    		find = true;

		i--;

	    }

	    i++;	    

    	}



    	// Assign min_gads to node and neighbours

    	if(find) {
	
    	    for( uint j = 0 ; j <= sp.second ; j++ ) {

    		if( (i+j) < Gads.size() )
    		    Gads[i+j] = make_pair( nodes[i+j], min_gads );

    		if( (i<j) > 0 )
    		    Gads[i-j] = make_pair( nodes[i-j], min_gads );	    
	    
    	    }

    	}
	

    }


    return Gads;

}
