#ifndef CARNAHANSTARLING_H
#define CARNAHANSTARLING_H

#include <EOS.H>

class CarnahanStarling : public EOS {


protected:

    
    /** Carnahan-Starling constants */

    scalar a;

    scalar b;


public:


    /** Constructor */

    CarnahanStarling( const std::string& dictName, const std::string& eqName );


    /** Default destructor */

    ~CarnahanStarling();


    /** EOS pressure */

    const scalar p_eos(const scalar& rho, const scalar& T) const;


    /** Pressure derivative */

    const scalar dp_dT(const scalar& rho, const scalar& T) const;    
    

};

#endif // CARNAHANSTARLING_H
