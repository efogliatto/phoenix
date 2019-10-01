#include <linearCavity.H>

using namespace std;


/** Constructor */
    
linearCavity::linearCavity( const std::string& bdname ) : cavityModel(bdname) {}


/** Destructor */

linearCavity::~linearCavity() {}


/** Compute adhesive coefficients */

const std::vector< std::pair<uint, scalar> > linearCavity::coeffs(const std::vector<uint>& nodes,
								  const std::vector< std::vector<uint> >& loc,
								  const std::vector< std::pair<uint,uint> >& spots) const {



}
