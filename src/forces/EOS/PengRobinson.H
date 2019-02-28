#ifndef PENGROBINSON_H
#define PENGROBINSON_H

#include <EOS.H>

class PengRobinson : public EOS {


protected:

    
    /** Peng-Robinson constants */

    scalar a;

    scalar b;

    scalar w;


public:


    /** Constructor */

    PengRobinson( const std::string& dictName, const std::string& eqName );


    /** Default destructor */

    ~PengRobinson();


    /** EOS pressure */

    const scalar p_eos(const scalar& rho, const scalar& T) const;


    /** Pressure derivative */

    const scalar dp_dT(const scalar& rho, const scalar& T) const;    
    

};

#endif // PENGROBINSON_H
