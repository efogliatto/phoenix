#ifndef GONGSRTEQ_H
#define GONGSRTEQ_H

#include <energyEquation.H>

#include <hsCreator.H>


class GongSRTEq : public energyEquation {


protected:

    
    /** Constants */

    scalar _kappa;

    scalar _Cv;

    

    /** Equation of state */

    EOS* eos;


    
    /** Thermal conductivity model */

    enum class thmodel { constCond, constDiff };

    thmodel thermal;

    

    /** Equilibrium in population space */

    const void eqPS( std::vector<scalar>& n, const uint& id ) const;
  



public:

    
    /** Default constructor */

    GongSRTEq( const std::string& name,
	     const latticeMesh& mesh_,
	     timeOptions& Time_,
	     pdfField& pdf_,
	     const scalarField& rho_,
	     const vectorField& U_,
	     scalarField& T_);


    /** Default destructor */

    ~GongSRTEq();



    /** Collision process */

    const void collision();


    /** Set pdf to equilibrium values */

    const void setEquilibrium();

    
    /** Update macroscopic temperature */

    const void updateMacroTemperature();


    /** Equilibrium in population space */

    const void eqPS( std::vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_, const scalar& hs ) const;    

    
    /** Thermal conductivity at node */

    const scalar thermalCond( const uint& id ) const;


    /** Diffusivity constant recovered at macroscopic level */

    const scalar diffusivityConstant(const uint id = 0) const;    
    
};


#endif // GONGSRTEQ_H
