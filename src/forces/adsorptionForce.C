#include <adsorptionForce.H>

using namespace std;


/** Constructor */

adsorptionForce::adsorptionForce( const string& dictName,
				  const string& eqName,
				  const latticeMesh& mesh )

    : _mesh(mesh) {


    // Create eos

    EOSCreator creator;

    eos = creator.create(dictName, eqName);
    
}


/** Destructor */

adsorptionForce::~adsorptionForce();
