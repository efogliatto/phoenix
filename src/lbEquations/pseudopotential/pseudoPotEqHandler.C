#include <pseudoPotEqHandler.H>

using namespace std;



/** Constructor */

pseudoPotEqHandler::pseudoPotEqHandler( const string& name,
					const latticeMesh& mesh_,
					timeOptions& Time_,
					pdfField& pdf_,
					scalarField& rho_,
					vectorField& U_,
					scalarField& T_ )

    : _mesh(mesh_),
      _rho(rho_) {


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


    // // First update density and forces

    // _equation->updateMacroDensity();

    // _equation->updateForces();


    // _equation->locateContactNodes();
    

    // Apply boundary conditions

    for(uint i = 0 ; i < _boundaries.size() ; i++) {

    	_boundaries[i]->update( _equation );

    }
    

}



/** Update macroscopic velocity */

const void pseudoPotEqHandler::updateMacroVelocity() {

    
    // Update forces

    _equation->updateForces();
    


    // Update forces at boundaries
    
    for(uint i = 0 ; i < _boundaries.size() ; i++) 
    	_boundaries[i]->updateIntForce( _equation );

    

    // Update macroscopic velocity
    
    _equation->updateMacroVelocity();

}



/** Update potential as scalar field  */

const void pseudoPotEqHandler::updatePotential( scalarField& phi ) {

    _equation->updatePotential(phi);

}




/** Compute and set pressure field */

const void pseudoPotEqHandler::pressure( const scalarField& phi, scalarField& p ) {

    _equation->pressure(phi, p);

}
