#ifndef EQFIXEDT_H
#define EQFIXEDT_H

#include <energyFixedT.H>

class eqFixedT : public energyFixedT {


public:


    /** Constructor */

    eqFixedT( const std::string& eqName,
	      const std::string& bdName,
	      const latticeMesh& mesh,
	      const scalarField& rho,
	      const scalarField& T,
	      const vectorField& U,
	      pdfField& pdf );


    /** Destructor */

    ~eqFixedT();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );

};

#endif // EQFIXEDT_H
