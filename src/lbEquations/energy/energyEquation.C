#include <energyEquation.H>

using namespace std;


/** Constructor */

energyEquation::energyEquation( const string& name,
				      const latticeMesh& mesh_,
				      timeOptions& Time_,
				      pdfField& pdf_,
				      const scalarField& rho_,
				      const vectorField& U_,
				      scalarField& T_)
    
    : lbEquation(name, mesh_, Time_, pdf_),
      rho(rho_),
      U(U_),
      T(T_) {
    


    // // Read Boundary conditions

    // dictionary dict("start/boundaries");

    // const map< string, vector<uint> >& bnd = mesh.boundaries();

    // ppBndCreator BndCreator;

    // for(map< string, vector<uint> >::const_iterator iter = bnd.begin(); iter != bnd.end(); ++iter)  {
	
    // 	string bdname = iter->first;

    // 	_boundaries.push_back( BndCreator.create(name, bdname, mesh.boundaryNodes(bdname) ) );

    // }

}


/** Default destructor */

energyEquation::~energyEquation() {}










/** Compute local density */

scalar energyEquation::localTemperature( const uint& id ) {

    scalar t(0);

    const uint q = mesh.lmodel()->q();
    
    for( uint k = 0 ; k < q ; k++ )
    	t += _pdf[id][k];


    return t;    

}





/** Collision process */

const void energyEquation::collision() {}



/** Update macroscopic density */

const void energyEquation::updateMacroTemperature() {

    for( uint i = 0 ; i < mesh.npoints() ; i++ )
	T[i] = localTemperature(i);

}






// /** Update boundaries */

// void energyEquation::updateBoundaries() {

//     for(uint i = 0 ; i < _boundaries.size() ; i++) {

// 	_boundaries[i]->update( mesh, _pdf, rho, T, U, &F );

//     }

// }