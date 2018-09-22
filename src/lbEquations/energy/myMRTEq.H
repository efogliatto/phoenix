#ifndef MYMRTEQ_H
#define MYMRTEQ_H

#include <energyEquation.H>

#include <hsCreator.H>


class myMRTEq : public energyEquation {


protected:

    
    /** Constants */

    scalar _a1;

    scalar _a2;


    /** Equilibrium in moment space */

    const void eqMS( std::vector<scalar>& m, const uint& id ) const;

    const void eqMS( std::vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_ ) const;    


    /** Equilibrium in population space */

    const void eqPS( std::vector<scalar>& n, const uint& id ) const;

    const void eqPS( std::vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_ ) const;
  



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


    /** Update boundary conditions */

    void updateBoundaries();
    
    
};


#endif // MYMRTEQ_H
