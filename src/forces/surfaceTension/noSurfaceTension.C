#include <noSurfaceTension.H>

using namespace std;


/** Constructor */

noSurfaceTension::noSurfaceTension( const string& dictName,
				    const string& eqName,
				    const latticeMesh& mesh )

    : surfaceTension(dictName, eqName, mesh) {}



/** Destructor */

noSurfaceTension::~noSurfaceTension() {}





/** Force at specific node */

const void noSurfaceTension::ST( const uint& i, const scalarField& rho, const scalarField& T, vector<scalar>& C, interactionForce* _fi, const std::vector<scalar>& Tau ) const {

    const uint q = _mesh.lmodel()->q();
    
    for(uint j = 0 ; j < q ; j++)
	C[j] = 0;

}
