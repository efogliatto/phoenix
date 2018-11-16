#ifndef FIXEDDENSITYBFORCE_H
#define FIXEDDENSITYBFORCE_H

#include <buoyantForce.H>

class fixedDensityBForce : public buoyantForce {

protected:

    /** Reference density */

    scalar _rhoRef;


public:

    
    /** Constructor */

    fixedDensityBForce( const std::string& dictName, const std::string& eqName );


    /** Destructor */

    ~fixedDensityBForce();


    /** Force */

    const std::vector<scalar> force(const scalar& rho);

    const scalar force(const scalar& rho, const uint& i) { return (rho-_rhoRef)*_gravity[i]; }


    /** Update coefficients */

    const void update(const scalarField& rho) {}

};

#endif // FIXEDDENSITYBFORCE_H
