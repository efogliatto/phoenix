#ifndef LINEARCAVITY_H
#define LINEARCAVITY_H

#include <cavityModel.H>

class linearCavity : public cavityModel {

public:

    /** Constructor */
    
    linearCavity( const std::string& bdname );


    /** Destructor */

    ~linearCavity();

    
    /** Compute adhesive coefficients */

    const std::vector< std::pair<uint, scalar> > coeffs(const std::vector<uint>& nodes,
							const std::vector< std::vector<uint> >& loc,
							const std::vector< std::pair<uint,uint> >& spots) const;

};

#endif // LINEARCAVITY_H
