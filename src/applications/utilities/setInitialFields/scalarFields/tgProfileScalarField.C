#include "tgProfileScalarField.H"

#include <math.h>

using namespace std;



void tgProfileScalarField( scalarField& field, const latticeMesh& mesh, const string& fname, const string& sptype ) {


    // Dictionary properties
    
    dictionary dict("start/initialFields");

    string auxName( fname + "/internalField/" + sptype + "/" );

    

    // Start point

    scalar start = dict.lookUp<scalar>( auxName + "startPoint" );  
    
    
    // End point

    scalar end = dict.lookUp<scalar>( auxName + "endPoint" );


    // Inside value

    scalar inval = dict.lookUp<scalar>( auxName + "inside" );


    // Outside value

    scalar outval = dict.lookUp<scalar>( auxName + "outside" );


    // Width

    scalar W = dict.lookUp<scalar>( auxName + "width" );    
    

    


    // Move over points and asign values


    for( uint i = 0 ; i < mesh.local() ; i++ ) {

    	vector<int> point = mesh.latticePoint(i);

    	scalar p0( 2.0 * (point[2]-start) / W ),
    	       p1( 2.0 * (point[2]-end) / W);

    	field[i] = outval + (inval-outval)/2.0 * (tanh(p0) - tanh(p1));

    }



    // Sync values across processes

    field.sync();
    



}
