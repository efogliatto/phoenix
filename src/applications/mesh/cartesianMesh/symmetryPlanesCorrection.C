#include <symmetryPlanesCorrection.H>

using namespace std;


const bool isOnPlane( const vector<uint>& point, const string& plane, const scalar offset );


void symmetryPlanesCorrection( vector< vector<int> >& nb,
			       const string& plane,
			       const scalar offset,
			       const vector< vector<uint> >& points,
			       const unordered_map< string, vector<uint> >& boundaries,
			       const latticeModel* lbmodel ) {


    // Lattice constants

    vector<uint> symIdx = lbmodel->symIdx(plane);

    const uint q = lbmodel->q();


    // Move over all boundary points and check if they are located over selected plane

    for( auto& bnd : boundaries ) {

	for( auto& id : bnd.second ) {

	    if(  isOnPlane( points[id], plane, offset )  ) {

		for( uint k = 0 ; k < q ; k++ ) {

		    if( nb[id][k] == -1 )
			nb[id][k] = nb[id][ symIdx[k] ];

		}

	    }	   	    

	}

    }
    
    
    

}




const bool isOnPlane( const vector<uint>& point, const string& plane, const scalar offset ) {

    bool onPlane(false);
    

    if( plane == "OXY" ) {

	if( point[2] == (uint)offset ) {

	    onPlane = true;
	    
	}

    }

    else {

	if( plane == "OXZ" ) {

	    if( point[1] == (uint)offset ) {

		onPlane = true;
	    
	    }

	}


	else {

	    if( plane == "OYZ" ) {

		if( point[0] == (uint)offset ) {

		    onPlane = true;
	    
		}

	    }

	    else {

		if( plane == "OX" ) {

		    if( point[1] == (uint)offset ) {

			onPlane = true;
	    
		    }

		}


		else {

		    if( plane == "OY" ) {

			if( point[0] == (uint)offset ) {

			    onPlane = true;
	    
			}

		    }

		}		

	    }	    

	}	

    }

    return onPlane;

}
