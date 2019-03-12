#ifndef NOADS_H
#define NOADS_H

#include <adhesiveForce.H>

class noAds : public adhesiveForce {

    
public:

    
    /** Constructor */

    noAds( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, const interactionForce* Fi, timeOptions& Time );


    /** Destructor */

    ~noAds();



    /** Update force field */

    void update( scalarField& rho, scalarField& T ) {}

    
};

#endif // NOADS_H
