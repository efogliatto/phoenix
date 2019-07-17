#ifndef NEEXTFIXEDTSPOT_H
#define NEEXTFIXEDTSPOT_H

#include <NEExtFixedT.H>

#include <sphere.H>


class NEExtFixedTSpot : public NEExtFixedT {


public:

    
    /** Constructor */

    NEExtFixedTSpot( const std::string& eqName,
		     const std::string& bdName,
		     const latticeMesh& mesh,
		     const scalarField& rho,
		     const scalarField& T,
		     const vectorField& U,
		     pdfField& pdf );


    /** Destructor */

    ~NEExtFixedTSpot();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );


};

#endif // NEEXTFIXEDTSPOT_H
