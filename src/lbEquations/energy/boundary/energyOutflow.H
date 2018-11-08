#ifndef ENERGYOUTFLOW_H
#define ENERGYOUTFLOW_H

#include <energyBndCond.H>

class energyOutflow : public energyBndCond {

protected:

    
    /** Pre-allocated neighbour node ids */

    std::vector<uint> _nbid;


    /** Boundary normal */

    std::vector< std::vector<int> > _normal;



public:


    /** Constructor */

    energyOutflow( const std::string& eqName,
		   const std::string& bdName,
		   const latticeMesh& mesh,
		   const scalarField& rho,
		   const scalarField& T,
		   const vectorField& U,
		   pdfField& pdf );


    /** Destructor */

    ~energyOutflow();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq );

};

#endif // ENERGYOUTFLOW_H
