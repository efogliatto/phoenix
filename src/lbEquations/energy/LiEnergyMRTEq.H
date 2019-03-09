#ifndef LIENERGYMRTEQ_H
#define LIENERGYMRTEQ_H

#include <energyEquation.H>

#include <hsCreator.H>


class LiEnergyMRTEq : public energyEquation {


protected:

    
    /** Constants */

    scalar _kappa;

    scalar _Cv;


    /** Equilibrium in moment space */

    const void eqMS( std::vector<scalar>& m, const uint& id ) const;

    const void eqMS( std::vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_ ) const;    


    /** Equilibrium in population space */

    const void eqPS( std::vector<scalar>& n, const uint& id ) const;
  



public:

    
    /** Default constructor */

    LiEnergyMRTEq( const std::string& name,
		   const latticeMesh& mesh_,
		   timeOptions& Time_,
		   pdfField& pdf_,
		   const scalarField& rho_,
		   const vectorField& U_,
		   scalarField& T_);


    /** Default destructor */

    ~LiEnergyMRTEq();



    /** Collision process */

    const void collision();


    /** Set pdf to equilibrium values */

    const void setEquilibrium();


    /** Compute local temperature */

    const scalar localTemperature( const uint& id ) const { return energyEquation::localTemperature(id); }

    
    /** Update macroscopic temperature */

    const void updateMacroTemperature();


    /** Equilibrium in population space */

    const void eqPS( std::vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_, const scalar& hs ) const;    

    
    /** Thermal conductivity at node */

    virtual const scalar thermalCond( const uint& id ) const;
    
};


#endif // LIENERGYMRTEQ_H
