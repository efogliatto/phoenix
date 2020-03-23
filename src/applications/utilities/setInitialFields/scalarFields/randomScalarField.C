#include "randomScalarField.H"

#include <random>

using namespace std;



void randomScalarField( scalarField& field, const latticeMesh& mesh, const string& fname, const string& sptype ) {


    dictionary dict("start/initialFields");

    string auxName( fname + "/internalField/" + sptype + "/" );

    
   

    // Inside Value

    scalar in = dict.lookUp<scalar>( auxName + "value" );

    scalar pert = dict.lookUpOrDefault<scalar>( auxName + "perturbation", 1 );




    // Random seeds
    
    uniform_real_distribution<scalar> unif( (100.0-pert)/100.0, (100.0+pert)/100.0);

    default_random_engine re;

    re.seed( mesh.local() );
    

    

    // Move over points and asign values

    for( uint i = 0 ; i < mesh.local() ; i++ )
	field[i] = unif(re) * in;




    // Sync values across processes

    field.sync();
    



}
