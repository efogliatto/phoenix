#ifndef EQFIXEDTSPOT_H
#define EQFIXEDTSPOT_H

#include <eqFixedT.H>

#include <sphere.H>


class eqFixedTSpot : public eqFixedT {


public:

    
    /** Constructor */

    eqFixedTSpot( const std::string& eqName,
		  const std::string& bdName,
		  const latticeMesh& mesh,
		  const scalarField& rho,
		  const scalarField& T,
		  const vectorField& U,
		  pdfField& pdf );


    /** Destructor */

    ~eqFixedTSpot();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );


};

#endif // EQFIXEDTSPOT_H
