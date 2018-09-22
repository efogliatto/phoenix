#include <energyEqHandler.H>

using namespace std;


/** Constructor */

energyEqHandler::energyEqHandler(const std::string& name,
				 const latticeMesh& mesh_,
				 timeOptions& Time_,
				 pdfField& pdf_,
				 const scalarField& rho_,
				 const vectorField& U_,
				 scalarField& T_) {

    
    // Energy MRT equation

    EEquation EEq;

    _equation = EEq.create(name, mesh_, Time_, pdf_, rho_, U_, T_);



    // Read Boundary conditions

    dictionary dict("start/boundaries");

    const map< string, vector<uint> >& bnd = mesh_.boundaries();

    energyBndCreator BndCreator;

    for(map< string, vector<uint> >::const_iterator iter = bnd.begin(); iter != bnd.end(); ++iter)  {
	
    	string bdname = iter->first;

    	_boundaries.push_back(   BndCreator.create(name, bdname, mesh_, rho_, T_, U_, pdf_)   );

    }
    

}



/** Destructor */

energyEqHandler::~energyEqHandler() {}


/** Update boundaries */

const void energyEqHandler::updateBoundaries() {

    for(uint i = 0 ; i < _boundaries.size() ; i++) {

    	_boundaries[i]->update( _equation );

    }

}
