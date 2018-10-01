#ifndef AVGDENSITYBFORCE_H
#define AVGDENSITYBFORCE_H

#include <buoyantForce.H>

class avgDensityBForce : public buoyantForce {

protected:

    /** Reference density */

    scalar _rhoRef;


public:

    
    /** Constructor */

    avgDensityBForce( const std::string& dictName, const std::string& eqName );


    /** Destructor */

    ~avgDensityBForce();


    /** Force */

    const std::vector<scalar> force(const scalar& rho);

    const scalar force(const scalar& rho, const uint& i) { return (rho-_rhoRef)*_gravity[i]; }


    /** Update coefficients */

    const void update(const scalarField& rho);

};

#endif // AVGDENSITYBFORCE_H
