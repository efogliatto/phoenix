#include <linearCavity.H>

using namespace std;


/** Constructor */
    
linearCavity::linearCavity( const std::string& bdname ) : cavityModel(bdname) {}


/** Destructor */

linearCavity::~linearCavity() {}


/** Compute adhesive coefficients */

const std::vector< std::pair<uint, scalar> > linearCavity::coeffs(const std::vector<uint>& nodes,
								  const std::vector< std::vector<uint> >& loc,
								  const std::vector< std::pair<uint,uint> >& spots) const {

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

	    scalar delta = (max_gads - min_gads) / sp.second;	    
	
    	    for( uint j = 0 ; j <= sp.second ; j++ ) {

    		if( (i+j) < Gads.size() )
    		    Gads[i+j] = make_pair( nodes[i+j], min_gads + delta*j );

    		if( (i-j) > 0 )
    		    Gads[i-j] = make_pair( nodes[i-j], min_gads + delta*j );	    
	    
    	    }

    	}
	

    }


    return Gads;

}
