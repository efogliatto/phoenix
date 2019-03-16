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


    _equation->updateContactNodes( pseudoPotEqHandler::contactLine() );
    

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





/** Check for contact points over surfaces */

vector<uint> pseudoPotEqHandler::contactLine() const {


    // Nodes with contact points
    
    vector<uint> contactNodes;


    // Move over boundaries
    
    for(uint i = 0 ; i < _boundaries.size() ; i++) {	
	    
	
	const vector<uint>& bnodes = _boundaries[i]->bdNodes();       

	vector<scalar> rhoDev( bnodes.size() );


	// Compute derivatives

	for( uint j = 0 ; j < bnodes.size() ; j++ ) {

	    if( j == 0 ) {

		rhoDev[j] = _rho.at(bnodes[j+1]) - _rho.at(bnodes[j]);

	    }

	    else {

		if( j == bnodes.size()-1 ) {

		    rhoDev[j] = _rho.at(bnodes[j]) - _rho.at(bnodes[j-1]);

		}

		else {

		    rhoDev[j] = 0.5 * ( _rho.at(bnodes[j+1]) - _rho.at(bnodes[j-1])  );

		}

	    }

	}


       

	// Assign contact point with maximum local derivative

	for( uint j = 1 ; j < bnodes.size() - 1 ; j++ ) {

	    if(  ( rhoDev[j] > rhoDev[j-1] )  &&  ( rhoDev[j] >= rhoDev[j+1] )  ) {

		if(j < _mesh.local())
		    contactNodes.push_back( bnodes[j] );

	    }

	}



    }




    return contactNodes;
    

}
