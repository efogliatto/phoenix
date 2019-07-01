#include "neq.H"

using namespace std;

scalar neq( const latticeMesh& mesh, const scalarField& T, const vectorField& U, const scalar alpha_1, const scalar alpha_2, const uint id, const uint vid ) {


    scalar a(0);

    switch(vid) {

    case 0:

	a = T.at(id);
	
	break;

	
    case 1:

	a = alpha_1 * T.at(id);
	
	break;


    case 2:

	a = alpha_2 * T.at(id);
	
	break;


    case 3:

	a = T.at(id) * U.at(id)[0];
	
	break;


    case 4:

	a = -T.at(id) * U.at(id)[0];
	
	break;


    case 5:

	a = T.at(id) * U.at(id)[1];
	
	break;


    case 6:

	a = -T.at(id) * U.at(id)[1];
	
	break;


    default:
	
	break;	

    }


    
    
    
    return a;
    
}
