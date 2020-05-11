#ifndef ENERGYNORMALHEATFLUXSPOT_H
#define ENERGYNORMALHEATFLUXSPOT_H

#include <energyNormalHeatFlux.H>

#include <sphere.H>


class energyNormalHeatFluxSpot : public energyNormalHeatFlux { 
    

public:

    
    /** Constructor */

    energyNormalHeatFluxSpot( const std::string& eqName,
			      const std::string& bdName,
			      const latticeMesh& mesh,
			      const scalarField& rho,
			      const scalarField& T,
			      const vectorField& U,
			      pdfField& pdf );


    /** Destructor */

    ~energyNormalHeatFluxSpot();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );


};

#endif // ENERGYNORMALHEATFLUXSPOT_H
