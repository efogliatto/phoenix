#ifndef LIHS_H
#define LIHS_H

#include <heatSource.H>


class liHS : public heatSource {


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

    liHS( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~liHS();
    

    /** Update source field */

    void update( const scalarField& rho, const scalarField& T, const vectorField& U );
    

};

#endif // LIHS_H
