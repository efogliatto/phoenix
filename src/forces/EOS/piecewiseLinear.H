#ifndef PIECEWISELINEAR_H
#define PIECEWISELINEAR_H

#include <EOS.H>

class piecewiseLinear : public EOS {


protected:

    
    /** Van Der Waals constants */

    scalar _rho_1;

    scalar _rho_2;


public:


    /** Constructor */

    piecewiseLinear( const std::string& dictName, const std::string& eqName );


    /** Default destructor */

    ~piecewiseLinear();


    /** EOS pressure */

    const scalar p_eos(const scalar& rho, const scalar& T) const;


    /** Pressure derivative */

    const scalar dp_dT(const scalar& rho, const scalar& T) const;    
    

};

#endif // PIECEWISELINEAR_H
