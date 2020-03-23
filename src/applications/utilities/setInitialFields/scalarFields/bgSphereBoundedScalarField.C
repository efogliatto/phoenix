#include "bgSphereBoundedScalarField.H"

#include <math.h>

using namespace std;



void bgSphereBoundedScalarField( scalarField& field, const latticeMesh& mesh, const string& fname, const string& sptype ) {


    dictionary dict("start/initialFields");

    string auxName( fname + "/internalField/" + sptype + "/" );

    

    // Centre

    vector<scalar> centre = dict.lookUp< vector<scalar> >( auxName + "centre" );	
    

    // Inside Value

    scalar in = dict.lookUp<scalar>( auxName + "inside" );


    // Radius

    scalar radius = dict.lookUp<scalar>( auxName + "radius" );

    
    // Width

    scalar width = dict.lookUpOrDefault<scalar>( auxName + "width", 1 );    

    



    // Move over points and asign values


    for( uint i = 0 ; i < mesh.local() ; i++ ) {

	vector<int> point = mesh.latticePoint(i);


	scalar r = ( point[0] - centre[0] )  *  ( point[0] - centre[0] ) +
	           ( point[1] - centre[1] )  *  ( point[1] - centre[1] ) + 
	           ( point[2] - centre[2] )  *  ( point[2] - centre[2] );	


	field[i] = 0.5*( in + field[i] )    +   0.5*( field[i] - in ) * tanh( (2*(sqrt(r) - radius)) / width )  ;


    }



    // Sync values across processes

    field.sync();
    



}
