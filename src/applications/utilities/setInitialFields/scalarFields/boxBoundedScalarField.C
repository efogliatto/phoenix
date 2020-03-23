#include "boxBoundedScalarField.H"

using namespace std;



void boxBoundedScalarField( scalarField& field, const latticeMesh& mesh, const string& fname, const string& sptype ) {


    dictionary dict("start/initialFields");

    string auxName( fname + "/internalField/" + sptype + "/" );

    

    // Min corner

    vector<scalar> min = dict.lookUp< vector<scalar> >( auxName + "minPoint" );
	

    // Max corner

    vector<scalar> max = dict.lookUp< vector<scalar> >( auxName + "maxPoint" );
    

    // Inside Value

    scalar in = dict.lookUp<scalar>( auxName + "inside" );


    // Outside Value

    scalar out = dict.lookUp<scalar>( auxName + "outside" );



    // Move over points and asign values


    for( uint i = 0 ; i < mesh.local() ; i++ ) {

	vector<int> point = mesh.latticePoint(i);
	

    	if( ( min[0] <= point[0] )   &&
    	    ( max[0] >= point[0] )   &&
    	    ( min[1] <= point[1] )   &&
    	    ( max[1] >= point[1] )   &&
    	    ( min[2] <= point[2] )   &&
    	    ( max[2] >= point[2] )  ) {
	
    	    field[i] = in;

    	}

    	else {

    	    field[i] = out;

    	}

    }



    // Sync values across processes

    field.sync();
    



}
