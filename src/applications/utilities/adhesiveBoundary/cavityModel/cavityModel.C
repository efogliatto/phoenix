#include <cavityModel.H>

using namespace std;


/** Constructor */

cavityModel::cavityModel( const std::string& bdname ) {

    
    // Read information

    dictionary dict("properties/adhesiveProperties");

    max_gads = dict.lookUpOrDefault<scalar>("Boundaries/" + bdname + "/CavityModel/max_gads", 0);

    min_gads = dict.lookUp<scalar>("Boundaries/" + bdname + "/CavityModel/min_gads");


}


/** Destructor */

cavityModel::~cavityModel() {}


/** Compute adhesive coefficients */

const std::vector< std::pair<uint, scalar> > cavityModel::coeffs(const std::vector<uint>& nodes,
								 const std::vector< std::vector<uint> >& loc,
								 const std::vector< std::pair<uint,uint> >& spots) const {

    return { std::make_pair(0,0) };

}
