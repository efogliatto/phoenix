#ifndef PPPERIODIC_H
#define PPPERIODIC_H

#include <ppBndCond.H>

class ppPeriodic : public ppBndCond {

    
public:


    /** Constructor */

    ppPeriodic( const std::vector<uint>& nodes );


    /** Destructor */

    ppPeriodic();


    /** Update pdf field */

    void update( const latticeMesh& mesh, pdfField& pdf, const scalarField& rho, const scalarField& T, const vectorField& U, const pseudoPotForce* F );

};

#endif // PPPERIODIC_H
