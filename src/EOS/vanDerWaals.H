#ifndef VANDERWAALS_H
#define VANDERWAALS_H

#include <EOS.H>

class vanDerWaals : public EOS {


protected:

    
    /** Van Der Waals constants */

    scalar a;

    scalar b;


public:


    /** Default constructor */

    vanDerWaals( const std::string fname );


    /** Default destructor */

    ~vanDerWaals();


    /** EOS pressure */

    const scalar p_eos(const scalar& rho, const scalar& T) const;
    

};

#endif // VANDERWAALS_H
