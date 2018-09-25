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


    /** Boundary type */

    const std::string _bdtype;    
   


public:

    /** Constructor */

    ppBndCond( const latticeMesh& mesh,
	       const scalarField& rho,
	       const scalarField& T,
	       const vectorField& U,
	       pdfField& pdf,
	       const std::string& bdName,
	       const std::string& bdType);


    /** Destructor */

    ~ppBndCond();


    /** Boundary type */

    const std::string& type() const { return _bdtype; }
    

    /** Update pdf field */

    virtual void update( const pseudoPotEquation* ppeq ) = 0;


    /** Update interaction force */

    virtual const void updateIntForce( pseudoPotEquation* ppeq ) const {}

};



#endif // PPBNDCOND_H
