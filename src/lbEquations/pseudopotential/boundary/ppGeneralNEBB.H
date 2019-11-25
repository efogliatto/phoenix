#ifndef PPGENERALNEBB_H
#define PPGENERALNEBB_H

#include <ppWetNodeBnd.H>

class ppGeneralNEBB : public ppWetNodeBnd {

protected:

    
    // Transformation coefficients

    std::map< latticeMesh::normalType, std::vector<uint> > trIndex;

    std::map< latticeMesh::normalType, std::vector<int> > prDirIndex;

    std::map< latticeMesh::normalType, std::vector<int> > prDirSign;        
    
    
public:


    /** Constructor */

    ppGeneralNEBB( const std::string& eqName,
		   const std::string& bdName,
		   const latticeMesh& mesh,
		   const scalarField& rho,
		   const scalarField& T,
		   const vectorField& U,
		   pdfField& pdf );


    /** Destructor */

    ~ppGeneralNEBB();


    /** Update pdf field */

    void update( const pseudoPotEquation* ppeq );

};

#endif // PPGENERALNEBB_H
