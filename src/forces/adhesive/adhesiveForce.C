#include <adhesiveForce.H>

using namespace std;


/** Constructor */

adhesiveForce::adhesiveForce( const string& dictName,
				  const string& eqName,
				  const latticeMesh& mesh )

    : _mesh(mesh) {


    // Create eos

    EOSCreator creator;

    eos = creator.create(dictName, eqName);
    
}


/** Destructor */

adhesiveForce::~adhesiveForce() {}
