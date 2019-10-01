#ifndef CAVITYMODEL_H
#define CAVITYMODEL_H

#include <dictionary.H>

class cavityModel {

protected:

    /** Maximun GAds */

    scalar max_gads;


    /** Minimun GAds */

    scalar min_gads;    
    
    
public:

    /** Constructor */

    cavityModel( const std::string& bdname );


    /** Destructor */

    ~cavityModel();


    /** Compute adhesive coefficients */

    virtual const std::vector< std::pair<uint, scalar> > coeffs(const std::vector<uint>& nodes,
								const std::vector< std::vector<uint> >& loc,
								const std::vector< std::pair<uint,uint> >& spots) const = 0;

};

#endif // CAVITYMODEL_H
