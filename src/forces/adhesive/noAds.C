#include <noAds.H>

using namespace std;


/** Constructor */

noAds::noAds( const string& dictName,
	      const string& eqName,
	      const latticeMesh& mesh,
	      const interactionForce* Fi,
	      timeOptions& Time )

    : adhesiveForce(dictName, eqName, mesh, Fi, Time) {}



/** Destructor */

noAds::~noAds() {}
