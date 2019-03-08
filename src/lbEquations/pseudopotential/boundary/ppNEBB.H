#ifndef PPNEBB_H
#define PPNEBB_H

#include <ppWetNodeBnd.H>

class ppNEBB : public ppWetNodeBnd {

public:


    /** Constructor */

    ppNEBB( const std::string& eqName,
	    const std::string& bdName,
	    const latticeMesh& mesh,
	    const scalarField& rho,
	    const scalarField& T,
	    const vectorField& U,
	    pdfField& pdf );


    /** Destructor */

    ~ppNEBB();


    /** Update pdf field */

    void update( const pseudoPotEquation* ppeq );

};

#endif // PPNEBB_H
