#ifndef PPOUTFLOW_H
#define PPOUTFLOW_H

#include <ppFixedU.H>

class ppOutflow : public ppFixedU {


protected:


    /** Pre-allocated neighbour node ids */

    std::vector<uint> _nbid;


public:

    
    /** Constructor */

    ppOutflow( const std::string& eqName,
	       const std::string& bdName,
	       const latticeMesh& mesh,
	       const scalarField& rho,
	       const scalarField& T,
	       const vectorField& U,
	       pdfField& pdf );


    /** Destructor */

    ~ppOutflow();


    /** Update pdf field */

    virtual void update( const pseudoPotEquation* ppeq );


};

#endif // PPOUTFLOW_H
