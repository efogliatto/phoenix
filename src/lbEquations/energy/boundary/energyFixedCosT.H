#ifndef ENERGYFIXEDCOST_H
#define ENERGYFIXEDCOST_H

#include <energyFixedT.H>

#include <sphere.H>


class energyFixedCosT : public energyFixedT {
  

public:

    
    /** Constructor */

    energyFixedCosT( const std::string& eqName,
		     const std::string& bdName,
		     const latticeMesh& mesh,
		     const scalarField& rho,
		     const scalarField& T,
		     const vectorField& U,
		     pdfField& pdf );


    /** Destructor */

    ~energyFixedCosT();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );


};

#endif // ENERGYFIXEDCOST_H
