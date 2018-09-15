#ifndef PPFIXEDU_H
#define PPFIXEDU_H

#include <ppBndCond.H>

class ppFixedU : public ppBndCond {

protected:

    /** Velocity value at boundary */

    std::vector< std::vector<scalar> > _bndVal;


public:


    /** Constructor */

    ppFixedU( const std::string& eqName, const std::string& bdName, const std::vector<uint>& nodes );


    /** Destructor */

    ppFixedU();


    /** Update pdf field */

    void update( const latticeMesh& mesh, pdfField& pdf, const scalarField& rho, const scalarField& T, const vectorField& U, const pseudoPotForce* F );

};

#endif // PPFIXEDU_H
