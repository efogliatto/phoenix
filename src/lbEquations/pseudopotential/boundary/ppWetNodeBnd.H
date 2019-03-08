#ifndef PPWETNODEBND_H
#define PPWETNODEBND_H

#include <ppBndCond.H>

class ppWetNodeBnd : public ppBndCond {

protected:

    /** Velocity value at boundary */

    std::vector< std::vector<scalar> > _bndVal;


    /** Pre-allocated neighbour node ids */

    std::vector<uint> _nbid;


    /** Boundary normal type */

    enum class normalType{ X0, X1, Y0, Y1, Z0, Z1, UNDEF };
    

    /** Boundary normal */

    std::vector<normalType> _normal;


public:


    /** Constructor */

    ppWetNodeBnd( const std::string& eqName,
		  const std::string& bdName,
		  const std::string& bdType,		  
		  const latticeMesh& mesh,
		  const scalarField& rho,
		  const scalarField& T,
		  const vectorField& U,
		  pdfField& pdf );


    /** Destructor */

    ~ppWetNodeBnd();


    /** Update pdf field */

    virtual void update( const pseudoPotEquation* ppeq ) = 0;


};

#endif // PPWETNODEBND_H
