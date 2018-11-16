#ifndef NOADS_H
#define NOADS_H

#include <adhesiveForce.H>

class noAds : public adhesiveForce {

    
public:

    
    /** Constructor */

    noAds( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~noAds();


    /** Force at specific node */

    const void force( const uint& i, const scalarField& rho, const scalarField& T, scalar f[3] ) const;

    
};

#endif // NOADS_H
