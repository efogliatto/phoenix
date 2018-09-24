#ifndef PPFIXEDU_H
#define PPFIXEDU_H

#include <ppBndCond.H>

class ppFixedU : public ppBndCond {

protected:

    /** Velocity value at boundary */

    std::vector< std::vector<scalar> > _bndVal;


public:


    /** Constructor */

    ppFixedU( const std::string& eqName,
	      const std::string& bdName,
	      const latticeMesh& mesh,
	      const scalarField& rho,
	      const scalarField& T,
	      const vectorField& U,
	      pdfField& pdf );


    /** Destructor */

    ~ppFixedU();


    /** Update pdf field */

    void update( const pseudoPotEquation* ppeq );

};

#endif // PPFIXEDU_H
