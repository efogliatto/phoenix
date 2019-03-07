#ifndef ENERGYFIXEDT_H
#define ENERGYFIXEDT_H

#include <energyBndCond.H>

class energyFixedT : public energyBndCond {

protected:

    /** Velocity value at boundary */

    std::vector<scalar> _bndVal;


    /** Perturbation (percentage) */

    scalar _pert;


    /** Pre-allocated neighbour node ids */

    std::vector<uint> _nbid;


    /** Boundary normal */

    std::vector< std::vector<int> > _normal;        



public:


    /** Constructor */

    energyFixedT( const std::string& eqName,
		  const std::string& bdName,
		  const latticeMesh& mesh,
		  const scalarField& rho,
		  const scalarField& T,
		  const vectorField& U,
		  pdfField& pdf );


    /** Destructor */

    ~energyFixedT();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );

};

#endif // ENERGYFIXEDT_H
