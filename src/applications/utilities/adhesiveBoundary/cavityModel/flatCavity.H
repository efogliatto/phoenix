#ifndef FLATCAVITY_H
#define FLATCAVITY_H

#include <cavityModel.H>

class flatCavity : public cavityModel {

public:

    /** Constructor */
    
    flatCavity( const std::string& bdname );


    /** Destructor */

    ~flatCavity();


    /** Compute adhesive coefficients */

    const std::vector< std::pair<uint, scalar> > coeffs(const std::vector<uint>& nodes,
							const std::vector< std::vector<uint> >& loc,
							const std::vector< std::pair<uint,uint> >& spots) const;


};

#endif // FLATCAVITY_H
