#ifndef RHOPIECEWISELINEAR_H
#define RHOPIECEWISELINEAR_H

#include <relaxModel.H>

#include <vector>


class rhoPieceWiseLinear : public relaxModel {

protected:

    
    /** Relaxation coefficients */
    
    std::vector<scalar> _Tau_v;

    std::vector<scalar> _Tau_l;


    /** Density limits */

    scalar _rho_v;

    scalar _rho_l;


public:

    /** Default constructor */

    rhoPieceWiseLinear( const std::string& entry );


    /** Destructor */

    ~rhoPieceWiseLinear();


    /** Density-dependent relaxation factor */

    const scalar tau( const scalar rho, const uint i );
    

};

#endif // RHOPIECEWISELINEAR_H
