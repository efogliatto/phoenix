#ifndef ENERGYEQUATION_H
#define ENERGYEQUATION_H

#include <lbEquation.H>

#include <hsCreator.H>



class energyEquation : public lbEquation {


protected:

    
    /** Reference to external density field */

    const scalarField& rho;


    /** Reference to external velocity field */

    const vectorField& U;


    /** Reference to external temperature field */

    scalarField& T;


    /** Heat sources for macro equation */

    heatSource* _hs;
  



public:

    
    /** Default constructor */

    energyEquation( const std::string& name,
		    const latticeMesh& mesh_,
		    timeOptions& Time_,
		    pdfField& pdf_,
		    const scalarField& rho_,
		    const vectorField& U_,
		    scalarField& T_);


    /** Default destructor */

    ~energyEquation();



    /** Collision process */

    virtual const void collision() = 0;


    /** Update macroscopic temperature */

    virtual const void updateMacroTemperature();


    /** Set pdf to equilibrium values */

    virtual const void setEquilibrium() = 0;


    /** Compute local temperature */

    virtual const scalar localTemperature( const uint& id ) const;

    
    /** Forced interface for equilibrium distribution */

    virtual const void eqPS( std::vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_ ) const = 0;    
    
    
};


#endif // ENERGYEQUATION_H
