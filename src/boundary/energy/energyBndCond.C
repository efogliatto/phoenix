#include <energyBndCond.H>


/** Constructor */

energyBndCond::energyBndCond( const std::vector<uint>& nodes,
			      const scalarField& rho,
			      const scalarField& T,
			      const vectorField& U,
			      pdfField& pdf )    
    : _nodes(nodes),
      _rho(rho),
      _T(T),
      _U(U),
      _pdf(pdf) {}



/** Destructor */

energyBndCond::~energyBndCond() {}
