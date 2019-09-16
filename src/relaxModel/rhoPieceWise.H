#ifndef RHOPIECEWISE_H
#define RHOPIECEWISE_H

#include <relaxModel.H>

#include <vector>


class rhoPieceWise : public relaxModel {

protected:

    
    /** Relaxation coefficients */
    
    std::vector<scalar> _Tau_v;

    std::vector<scalar> _Tau_l;


    /** Density limits */

    scalar _rho_lim;


public:

    /** Default constructor */

    rhoPieceWise( const std::string& entry );


    /** Destructor */

    ~rhoPieceWise();


    /** Density-dependent relaxation factor */

    const scalar tau( const scalar rho, const uint i );
    

};

#endif // RHOPIECEWISE_H
