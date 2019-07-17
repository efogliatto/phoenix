#ifndef NEEXTFIXEDT_H
#define NEEXTFIXEDT_H

#include <energyFixedT.H>

class NEExtFixedT : public energyFixedT {


public:


    /** Constructor */

    NEExtFixedT( const std::string& eqName,
		 const std::string& bdName,
		 const latticeMesh& mesh,
		 const scalarField& rho,
		 const scalarField& T,
		 const vectorField& U,
		 pdfField& pdf );


    /** Destructor */

    ~NEExtFixedT();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );

};

#endif // NEEXTFIXEDT_H
