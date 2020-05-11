#ifndef ENERGYNORMALHEATFLUX_H
#define ENERGYNORMALHEATFLUX_H

#include <energyFixedGradT.H>

class energyNormalHeatFlux : public energyFixedGradT {


protected:


    /** Temperature limit */

    scalar _tempLimit;
    
    

public:

    
    /** Constructor */

    energyNormalHeatFlux( const std::string& eqName,
			  const std::string& bdName,
			  const latticeMesh& mesh,
			  const scalarField& rho,
			  const scalarField& T,
			  const vectorField& U,
			  pdfField& pdf );


    /** Destructor */

    ~energyNormalHeatFlux();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );


};

#endif // ENERGYNORMALHEATFLUX_H
