#ifndef ENERGYBNDCOND_H
#define ENERGYBNDCOND_H

#include <EEquation.H>


class energyBndCond {

protected: 

    
    /** Reference to external fields */

    const std::vector<uint>& _nodes;

    const latticeMesh& _mesh;

    const scalarField& _rho;

    const scalarField& _T;

    const vectorField& _U;

    pdfField& _pdf;
    
      

public:

    /** Constructor */

    energyBndCond( const latticeMesh& mesh,
		   const scalarField& rho,
		   const scalarField& T,
		   const vectorField& U,
		   pdfField& pdf,
		   const std::string& bdName);


    /** Destructor */

    ~energyBndCond();


    /** Update pdf field */

    virtual void update( const energyEquation* eeq ) = 0;
    

};



#endif // ENERGYBNDCOND_H
