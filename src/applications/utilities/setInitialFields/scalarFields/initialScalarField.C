#include <initialScalarField.H>

#include "boxBoundedScalarField.H"

#include "uniformScalarField.H"

#include "randomScalarField.H"

#include "bgSphereBoundedScalarField.H"


using namespace std;


// Default constructor

initialScalarField::initialScalarField() {

    // Initialize types

    _spMapType["box"]        = ishape::box;
    _spMapType["random"]     = ishape::random;
    _spMapType["bgSphere"]   = ishape::bgsphere;
    _spMapType["uniform"]    = ishape::uniform;        

}



// Update field

void initialScalarField::updateField( scalarField& field, const latticeMesh& mesh, const string& fname, const string& sptype ) {


    // Read specific shape type from dictionary

    dictionary idict( "start/initialFields" );
    
    string sp = idict.lookUp<string>( fname + "/internalField/" + sptype + "/type" );


    
    
    if( _spMapType.find(sp) != _spMapType.end() ) {
	

    	switch( _spMapType[sp] ) {

    	case ishape::box:

	    boxBoundedScalarField(field, mesh, fname, sptype);

    	    break;

	    
    	case ishape::random:

	    randomScalarField(field, mesh, fname, sptype);
	    
    	    break;


    	case ishape::uniform:

	    uniformScalarField(field, mesh, fname, sptype);
	    
    	    break;	    


    	case ishape::bgsphere:

	    bgSphereBoundedScalarField(field, mesh, fname, sptype);
	    
    	    break;	    
	    
	    
    	}

    }


    else {
    
    	cout << endl << " [ERROR]  Initial shape " << sp << " not available" << endl << endl;

    	exit(1);

    }

    

}
