#ifndef ENERGYFIXEDT_H
#define ENERGYFIXEDT_H

#include <energyBndCond.H>

class energyFixedT : public energyBndCond {

protected:

    /** Velocity value at boundary */

    std::vector<scalar> _bndVal;


public:


    /** Constructor */

    energyFixedT( const std::string& eqName,
		  const std::string& bdName,
		  const std::vector<uint>& nodes,
		  const scalarField& rho,
		  const scalarField& T,
		  const vectorField& U,
		  pdfField& pdf );


    /** Destructor */

    ~energyFixedT();


    /** Update pdf field */

    virtual void update( const latticeMesh& mesh, std::function<void(std::vector<scalar>&, const scalar&, const std::vector<scalar>&)> eq );

};

#endif // ENERGYFIXEDT_H
