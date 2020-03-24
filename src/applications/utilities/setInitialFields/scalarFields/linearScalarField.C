#include "linearScalarField.H"

#include <math.h>

using namespace std;



void linearScalarField( scalarField& field, const latticeMesh& mesh, const string& fname, const string& sptype ) {


    // Dictionary properties
    
    dictionary dict("start/initialFields");

    string auxName( fname + "/internalField/" + sptype + "/" );

    

    // Start point

    vector<scalar> start = dict.lookUp< vector<scalar> >( auxName + "startPoint" );

    scalar stval = dict.lookUp<scalar>( auxName + "startValue" );
    
    
    
    // End point

    vector<scalar> end = dict.lookUp< vector<scalar> >( auxName + "endPoint" );

    scalar endval = dict.lookUp<scalar>( auxName + "endValue" );    
    


    // Compute linear constants

    scalar dmag(   (end[0]-start[0])*(end[0]-start[0])
	         + (end[1]-start[1])*(end[1]-start[1])
	         + (end[2]-start[2])*(end[2]-start[2])  );

    dmag = dmag*(endval-stval);

    scalar A[3] = { (end[0]-start[0]) / dmag,
                    (end[1]-start[1]) / dmag,
                    (end[2]-start[2]) / dmag };

    


    // Move over points and asign values


    for( uint i = 0 ; i < mesh.local() ; i++ ) {

	vector<int> point = mesh.latticePoint(i);

	scalar aux(0);

	for(uint j = 0 ; j < 3 ; j++)
	    aux += A[j] * (point[j] - start[j]) ;


	field[i] = aux + stval;


    }



    // Sync values across processes

    field.sync();
    



}
