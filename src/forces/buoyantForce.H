#ifndef BUOYANTFORCE_H
#define BUOYANTFORCE_H

#include <scalarField.H>

#include <vectorField.H>


class buoyantForce {

protected:

    /** Gravity magnitude */
    
    std::vector<scalar> _gravity;


public:

    /** Constructor */

    buoyantForce( const std::string& dictName, const std::string& eqName );


    /** Destructor */

    ~buoyantForce();


    /** Force */

    virtual const std::vector<scalar> force(const scalar& rho) = 0; 

};

#endif // BUOYANTFORCE_H
