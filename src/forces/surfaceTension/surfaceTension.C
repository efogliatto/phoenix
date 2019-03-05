#include <surfaceTension.H>

using namespace std;


/** Constructor */

surfaceTension::surfaceTension( const string& dictName,
				const string& eqName,
				const latticeMesh& mesh )

    : _mesh(mesh) {


    // // Create eos

    // EOSCreator creator;

    // eos = creator.create(dictName, eqName);
    
}


/** Destructor */

surfaceTension::~surfaceTension() {}
