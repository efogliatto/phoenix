#ifndef PPOUTFLOWWITHNEBB_H
#define PPOUTFLOWWITHNEBB_H

#include <ppGeneralNEBB.H>

class ppOutflowWithNEBB : public ppGeneralNEBB {


protected:


    // /** Pre-allocated neighbour node ids */

    // std::vector<uint> _nbid;


    // /** Boundary normal */

    // std::vector< std::vector<int> > _normal;

    
    // /** Old values of pdf */

    // std::vector< std::vector<scalar> > _oldPdf;


    // Lattice index for normal neighbours

    uint ln;


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
