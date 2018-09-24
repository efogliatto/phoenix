#ifndef PPPERIODIC_H
#define PPPERIODIC_H

#include <ppBndCond.H>

class ppPeriodic : public ppBndCond {

    
public:


    /** Constructor */

    ppPeriodic( const std::string& eqName,
		const std::string& bdName,
		const latticeMesh& mesh,
		const scalarField& rho,
		const scalarField& T,
		const vectorField& U,
		pdfField& pdf );


    /** Destructor */

    ~ppPeriodic();


    /** Update pdf field */

    void update( const pseudoPotEquation* ppeq );

};

#endif // PPPERIODIC_H
