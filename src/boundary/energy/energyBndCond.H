#ifndef ENERGYBNDCOND_H
#define ENERGYBNDCOND_H

#include <pdfField.H>

#include <scalarField.H>

#include <vectorField.H>

#include <functional>


class energyBndCond {

protected:

    
    /** Reference to node ids */

    const std::vector<uint> _nodes;

    
    /** Reference to external fields */

    const scalarField& _rho;

    const scalarField& _T;

    const vectorField& _U;

    pdfField& _pdf;
    
      

public:

    /** Constructor */

    energyBndCond( const std::vector<uint>& nodes,
		   const scalarField& rho,
		   const scalarField& T,
		   const vectorField& U,
		   pdfField& pdf);


    /** Destructor */

    ~energyBndCond();


    /** Update pdf field */

    virtual void update( const latticeMesh& mesh ) = 0;

};



#endif // ENERGYBNDCOND_H
