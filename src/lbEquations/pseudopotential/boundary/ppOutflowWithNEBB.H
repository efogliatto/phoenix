#ifndef PPOUTFLOWWITHNEBB_H
#define PPOUTFLOWWITHNEBB_H

#include <ppGeneralNEBB.H>

class ppOutflowWithNEBB : public ppGeneralNEBB {


protected:


    // Lattice index for normal neighbours

    uint ln;


    // Boundary mean normal

    std::vector<scalar> _normal;


    // Pre allocated density map

    std::map< uint, scalar > _rhow;


public:

    
    /** Constructor */

    ppOutflowWithNEBB( const std::string& eqName,
		       const std::string& bdName,
		       const latticeMesh& mesh,
		       const scalarField& rho,
		       const scalarField& T,
		       const vectorField& U,
		       pdfField& pdf);


    /** Destructor */

    ~ppOutflowWithNEBB();


    /** Update pdf field */

    virtual void update( const pseudoPotEquation* ppeq );


    // /** Update interaction force */
    
    // const void updateIntForce( pseudoPotEquation* ppeq ) const;



};

#endif // PPOUTFLOWWITHNEBB_H
