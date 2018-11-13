#ifndef ENERGYFIXEDTSPOT_H
#define ENERGYFIXEDTSPOT_H

#include <energyFixedT.H>

#include <sphere.H>


class energyFixedTSpot : public energyFixedT {


public:

    
    /** Constructor */

    energyFixedTSpot( const std::string& eqName,
		      const std::string& bdName,
		      const latticeMesh& mesh,
		      const scalarField& rho,
		      const scalarField& T,
		      const vectorField& U,
		      pdfField& pdf );


    /** Destructor */

    ~energyFixedTSpot();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );


};

#endif // ENERGYFIXEDTSPOT_H
