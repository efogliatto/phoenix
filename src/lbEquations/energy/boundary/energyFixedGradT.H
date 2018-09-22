#ifndef ENERGYFIXEDGRADT_H
#define ENERGYFIXEDGRADT_H

#include <energyFixedT.H>

class energyFixedGradT : public energyFixedT {


protected:

    /** Gradient value */

    scalar _grad;


    /** Pre-allocated neighbour node ids */

    std::vector<uint> _nbid;


public:

    
    /** Constructor */

    energyFixedGradT( const std::string& eqName,
		      const std::string& bdName,
		      const latticeMesh& mesh,
		      const scalarField& rho,
		      const scalarField& T,
		      const vectorField& U,
		      pdfField& pdf );


    /** Destructor */

    ~energyFixedGradT();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );


};

#endif // ENERGYFIXEDGRADT_H
