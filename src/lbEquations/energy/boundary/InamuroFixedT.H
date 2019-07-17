#ifndef INAMUROFIXEDT_H
#define INAMUROFIXEDT_H

#include <energyFixedT.H>

class InamuroFixedT : public energyFixedT {


public:


    /** Constructor */

    InamuroFixedT( const std::string& eqName,
		   const std::string& bdName,
		   const latticeMesh& mesh,
		   const scalarField& rho,
		   const scalarField& T,
		   const vectorField& U,
		   pdfField& pdf );


    /** Destructor */

    ~InamuroFixedT();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );

};

#endif // INAMUROFIXEDT_H
