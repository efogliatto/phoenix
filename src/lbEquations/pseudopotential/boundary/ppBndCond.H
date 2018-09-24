#ifndef PPBNDCOND_H
#define PPBNDCOND_H

#include <PPEquation.H>


class ppBndCond {

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

    ppBndCond( const latticeMesh& mesh,
	       const scalarField& rho,
	       const scalarField& T,
	       const vectorField& U,
	       pdfField& pdf,
	       const std::string& bdName);


    /** Destructor */

    ~ppBndCond();


    /** Update pdf field */

    virtual void update( const pseudoPotEquation* ppeq ) = 0;

};



#endif // PPBNDCOND_H
