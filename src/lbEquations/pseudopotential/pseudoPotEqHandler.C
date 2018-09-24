#include <pseudoPotEqHandler.H>

using namespace std;



/** Constructor */

pseudoPotEqHandler::pseudoPotEqHandler( const string& name,
					const latticeMesh& mesh_,
					timeOptions& Time_,
					pdfField& pdf_,
					scalarField& rho_,
					vectorField& U_,
					scalarField& T_ ) {


    // Energy MRT equation

    PPEquation PPEq;

    _equation = PPEq.create(name, mesh_, Time_, pdf_, rho_, U_, T_);



    // Read Boundary conditions

    dictionary dict("start/boundaries");

    const map< string, vector<uint> >& bnd = mesh_.boundaries();

    ppBndCreator BndCreator;

    for(map< string, vector<uint> >::const_iterator iter = bnd.begin(); iter != bnd.end(); ++iter)  {
	
    	string bdname = iter->first;

    	_boundaries.push_back(   BndCreator.create(name, bdname, mesh_, rho_, T_, U_, pdf_)   );

    }
    

}


/** Destructor */

pseudoPotEqHandler::~pseudoPotEqHandler() {}



/** Update boundaries */

const void pseudoPotEqHandler::updateBoundaries() {

    for(uint i = 0 ; i < _boundaries.size() ; i++) {

    	_boundaries[i]->update( _equation );

    }

}


