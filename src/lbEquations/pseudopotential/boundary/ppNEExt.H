#ifndef PPNEEXT_H
#define PPNEEXT_H

#include <ppWetNodeBnd.H>

class ppNEExt : public ppWetNodeBnd {

public:


    /** Constructor */

    ppNEExt( const std::string& eqName,
	     const std::string& bdName,
	     const latticeMesh& mesh,
	     const scalarField& rho,
	     const scalarField& T,
	     const vectorField& U,
	     pdfField& pdf );


    /** Destructor */

    ~ppNEExt();


    /** Update pdf field */

    void update( const pseudoPotEquation* ppeq );

};

#endif // PPNEEXT_H
