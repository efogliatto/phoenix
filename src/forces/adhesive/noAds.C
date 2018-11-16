#include <noAds.H>

using namespace std;


/** Constructor */

noAds::noAds( const string& dictName,
	      const string& eqName,
	      const latticeMesh& mesh )

    : adhesiveForce(dictName, eqName, mesh) {}



/** Destructor */

noAds::~noAds() {}


/** Force at specific node */

const void noAds::force( const uint& i, const scalarField& rho, const scalarField& T, scalar f[3] ) const {

    f[0] = 0;
    f[1] = 0;
    f[2] = 0;

}
