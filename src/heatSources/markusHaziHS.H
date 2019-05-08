#ifndef MARKUSHAZIHS_H
#define MARKUSHAZIHS_H

#include <heatSource.H>


class markusHaziHS : public heatSource {


protected:

    
    /** Equation of state */

    EOS* eos;

    
    /** Copy of relaxation coefficients */

    std::vector<scalar> _Tau;


    /** Constant specific heat */

    scalar _Cv;


    /** Constants */

    scalar _a1;

    scalar _a2;
    


public:


    /** Constructor */

    markusHaziHS( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~markusHaziHS();
    

    /** Update source field */

    void update( const scalarField& rho, const scalarField& T, const vectorField& U );


    /** Update source field with external temperature gradient */
    
    void update( const scalarField& rho, const scalarField& T, const vectorField& U, const vectorField& Tgrad );    
    

};

#endif // MARKUSHAZIHS_H
