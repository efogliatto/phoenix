#ifndef UNIFORMTAU_H
#define UNIFORMTAU_H

#include <relaxModel.H>

#include <vector>


class uniformTau : public relaxModel {

protected:

    
    /** Relaxation coefficients */
    
    std::vector<scalar> _Tau;


public:

    /** Default constructor */

    uniformTau( const std::string& entry );


    /** Destructor */

    ~uniformTau();


    /** Density-dependent relaxation factor */

    const scalar tau( const scalar rho, const uint i );
    

};

#endif // UNIFORMTAU_H
