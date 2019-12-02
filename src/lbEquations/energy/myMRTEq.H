#ifndef MYMRTEQ_H
#define MYMRTEQ_H

#include <energyEquation.H>

#include <hsCreator.H>


class myMRTEq : public energyEquation {


protected:


    /** Reference to temperature gradient field */

    vectorField Tgrad;


    /** Flag to enable temperature gradient calculation with local population */

    bool localGradT;
    
    
    /** Constants */

    scalar _a1;

    scalar _a2;

    scalar _Cv;


    /** Equilibrium in moment space */

    const void eqMS( std::vector<scalar>& m, const uint& id ) const;

    const void eqMS( std::vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_ ) const;    


    /** Equilibrium in population space */

    const void eqPS( std::vector<scalar>& n, const uint& id ) const;
  



public:

    
    /** Default constructor */

    myMRTEq( const std::string& name,
	     const latticeMesh& mesh_,
	     timeOptions& Time_,
	     pdfField& pdf_,
	     const scalarField& rho_,
	     const vectorField& U_,
	     scalarField& T_);


    /** Default destructor */

    ~myMRTEq();



    /** Collision process */

    const void collision();


    /** Set pdf to equilibrium values */

    const void setEquilibrium();


    /** Compute local temperature */

    const scalar localTemperature( const uint& id ) const { return energyEquation::localTemperature(id) + 0.5*_hs->source(id); }

    
    /** Update macroscopic temperature */

    const void updateMacroTemperature();


    /** Equilibrium in population space */

    const void eqPS( std::vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_, const scalar& hs ) const;    

    
    /** Thermal conductivity at node */

    virtual const scalar thermalCond( const uint& id ) const;


    /** Diffusivity constant recovered at macroscopic level */

    const scalar diffusivityConstant( const uint id = 0 ) const;    
    
};


#endif // MYMRTEQ_H
