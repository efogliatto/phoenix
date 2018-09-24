#include <ppBndCond.H>

using namespace std;


/** Constructor */

ppBndCond::ppBndCond( const latticeMesh& mesh,
		      const scalarField& rho,
		      const scalarField& T,
		      const vectorField& U,
		      pdfField& pdf,
		      const std::string& bdName )
    
    : _nodes( mesh.boundaryNodes(bdName) ),
      _mesh(mesh),
      _rho(rho),
      _T(T),
      _U(U),
      _pdf(pdf) {}


/** Destructor */

ppBndCond::~ppBndCond() {}
