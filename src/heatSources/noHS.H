#ifndef NOHS_H
#define NOHS_H

#include <heatSource.H>


class noHS : public heatSource {


protected:

    
    /** Equation of state */

    EOS* eos;

    
    /** Copy of relaxation coefficients */

    std::vector<scalar> _Tau;


    /** Constant specific heat */

    scalar _Cv;


    /** Constants */

    scalar _kappa;


    /** Thermal conductivity model */

    enum class thmodel { constCond, constDiff };

    thmodel thermal;
    


public:


    /** Constructor */

    noHS( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~noHS();
    

    /** Update source field */

    void update( const scalarField& rho, const scalarField& T, const vectorField& U );
    

};

#endif // NOHS_H
