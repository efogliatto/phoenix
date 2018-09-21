#ifndef EOS_H
#define EOS_H

#include <dataTypes.H>

#include <string>

#include <cmath>


class EOS {

    
protected:

    
    /** Equation of state name */

    const std::string name;   
    
    


public:


    /** Default constructor */

    EOS( const std::string& name_ );


    /** Default destructor */

    ~EOS();


    /** EOS pressure */

    virtual const scalar p_eos(const scalar& rho, const scalar& T) const = 0;

    
    /** Pressure derivative */

    virtual const scalar dp_dT(const scalar& rho, const scalar& T) const = 0;
    

};

#endif // EOS_H
