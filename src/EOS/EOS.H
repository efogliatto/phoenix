#ifndef EOS_H
#define EOS_H

#include <dataTypes.H>

#include <string>

#include <cmath>


class EOS {

    
protected:

    
    /** Equation of state name */

    const std::string name;

    
    /** Interaction strength */

    scalar _G;

    
    


public:


    /** Default constructor */

    EOS( const std::string& name_ );


    /** Default destructor */

    ~EOS();


    /** Pseudo-potential */

    const scalar potential( const scalar& p, const scalar& cs2, const scalar& rho, const scalar& T ) const {

	return sqrt(   2 * ( p - rho * cs2 ) / _G   );

    }


    /** EOS pressure */

    virtual const scalar p_eos(const scalar& rho, const scalar& T) const = 0;


    /** Interaction strength */

    const scalar& G() const { return _G; }

};

#endif // EOS_H
