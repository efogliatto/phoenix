#ifndef ENERGYPERIODIC_H
#define ENERGYPERIODIC_H

#include <energyBndCond.H>

class energyPeriodic : public energyBndCond {

protected:

    /** Velocity value at boundary */

    std::vector<scalar> _bndVal;


public:


    /** Constructor */

    energyPeriodic( const std::string& eqName,
		    const std::string& bdName,
		    const latticeMesh& mesh,		    
		    const scalarField& rho,
		    const scalarField& T,
		    const vectorField& U,
		    pdfField& pdf );


    /** Destructor */

    ~energyPeriodic();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );

};

#endif // ENERGYPERIODIC_H
