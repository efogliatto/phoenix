#ifndef PPBNDCOND_H
#define PPBNDCOND_H

#include <pdfField.H>

#include <pseudoPotForce.H>


class ppBndCond {

protected:

    
    /** Reference to node ids */

    const std::vector<uint> _nodes;


public:

    /** Constructor */

    ppBndCond( const std::vector<uint>& nodes );


    /** Destructor */

    ppBndCond();


    /** Update pdf field */

    virtual void update( const latticeMesh& mesh, pdfField& pdf, const scalarField& rho, const scalarField& T, const vectorField& U, const pseudoPotForce* F ) = 0;

};



#endif // PPBNDCOND_H
