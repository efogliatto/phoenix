#include "uniformScalarField.H"

using namespace std;



void uniformScalarField( scalarField& field, const latticeMesh& mesh, const string& fname, const string& sptype ) {


    dictionary dict("start/initialFields");

    string auxName( fname + "/internalField/" + sptype + "/" );

    
   

    // Inside Value

    scalar in = dict.lookUp<scalar>( auxName + "value" );



    // Move over points and asign values

    for( uint i = 0 ; i < mesh.local() ; i++ )
	field[i] = in;




    // Sync values across processes

    field.sync();
    



}
