#include <pseudoPotEquation.H>

#include <algorithm>

using namespace std;


/** Constructor */

pseudoPotEquation::pseudoPotEquation( const string& name,
				      const latticeMesh& mesh_,
				      timeOptions& Time_,
				      pdfField& pdf_,
				      scalarField& rho_,
				      vectorField& U_,
				      scalarField& T_)
    
    : lbEquation(name, mesh_, Time_, pdf_),
      rho(rho_),
      U(U_),
      T(T_),
      F("properties/macroProperties", name, mesh_, Time_, rho_, T_) {


    
    // Update Forces

    F.update(rho,T);


    // Read contact angles for each boundary

    dictionary dict("properties/macroProperties");

    string contact = dict.lookUpOrDefault<string>( name + "/ContactAngle/type", "none");

    if( contact == "geometric" ) {


	// Density on interface

	_rhoAvgInt = dict.lookUp<scalar>( name + "/ContactAngle/rhoAvg");

	

	// Read angles for each boundary

	const map< string, vector<uint> >& boundary = mesh.boundaries();

	map<string, scalar> angle;

	for( const auto& bd : boundary) {

	    scalar ang = dict.lookUpOrDefault<scalar>( name + "/ContactAngle/Theta/" + bd.first, 90);
	    
	    angle[bd.first] = ang * M_PI / 180.0;

	}



	// Create map for each boundary node

	for( const auto& bd : angle) {

	    for( const auto& id : boundary.at(bd.first) ) {

		_contactAngle[id] = bd.second;

	    }

	}
	

    }
    

}


/** Default destructor */

pseudoPotEquation::~pseudoPotEquation() {}





/** Compute local density */

const scalar pseudoPotEquation::localDensity( const uint& id ) const {

    scalar r(0);


    
    const uint q = mesh.lmodel()->q();
    
    for( uint k = 0 ; k < q ; k++ )
    	r += _pdf[id][k];



    

    // // Regular node

    // if( std::find(_contactNodes.begin(), _contactNodes.end(), id) == _contactNodes.end() ) {

    // 	const uint q = mesh.lmodel()->q();
    
    // 	for( uint k = 0 ; k < q ; k++ )
    // 	    r += _pdf[id][k];

    // }

    // else {


    // 	const uint q = mesh.lmodel()->q();
    
    // 	for( uint k = 0 ; k < q ; k++ )
    // 	    r += _pdf[id][k];
	

    // 	const vector< vector<int> >& nb = mesh.nbArray();

    // 	int neigh = nb[ nb[id][4] ][4];

    // 	if(neigh == -1)
    // 	    neigh = nb[id][4];

    // 	scalar delta = rho.at(nb[id][7]) - rho.at(nb[id][8]);

    // 	delta = abs(delta);       

    // 	r = rho.at(neigh) + tan( M_PI/2 - _contactAngle.at(id) ) * delta;

    // }
    

    
    
    return r;    

}



/** Compute local density with unknowns */

const scalar pseudoPotEquation::localDensityWithUnknowns( const uint& id, latticeMesh::normalType& ntype ) const {


    scalar rw(0);

    scalar Ft[3] = {0,0,0};

    
    // Hand coded

    switch( mesh.lmodel()->type() ) {

    case latticeModel::latticeType::D2Q9:


	// Total force at node
	
	F.total(Ft, id);
	

	switch(ntype) {

	case latticeMesh::normalType::Y1:
	    
	    rw = _pdf[id][0] + _pdf[id][1] + _pdf[id][3] + 2*(_pdf[id][2] + _pdf[id][5] + _pdf[id][6]) + 0.5 * Ft[1];

	    rw = rw / (1 + U.at(id)[1]);
	    
	    break;


	case latticeMesh::normalType::Y0:
	    
	    rw = _pdf[id][0] + _pdf[id][1] + _pdf[id][3] + 2*(_pdf[id][4] + _pdf[id][7] + _pdf[id][8]) - 0.5 * Ft[1];

	    rw = rw / (1 - U.at(id)[1]);
	    
	    break;


	case latticeMesh::normalType::X0:
	    
	    rw = _pdf[id][0] + _pdf[id][2] + _pdf[id][4] + 2*(_pdf[id][3] + _pdf[id][6] + _pdf[id][7]) - 0.5 * Ft[0];

	    rw = rw / (1 - U.at(id)[0]);
	    
	    break;


	case latticeMesh::normalType::X1:
	    
	    rw = _pdf[id][0] + _pdf[id][2] + _pdf[id][4] + 2*(_pdf[id][1] + _pdf[id][5] + _pdf[id][8]) + 0.5 * Ft[0];

	    rw = rw / (1 + U.at(id)[0]);
	    
	    break;	    
	    

	default:

	    rw = rho.at(id);
	    
	    break;

	}

	
	break;


    default:

	cout << "Local density with unknowns not yet implemented for this lattice type" << endl;

	exit(1);

	break;
	
    }

    

    return rw;
    

}




/** Compute local velocity */

const void pseudoPotEquation::localVelocity( vector<scalar>& v, const uint& id, const bool updDens ) const {
    

    // Local velocity
    
    vector<scalar> lv = {0,0,0};


    // Local density

    scalar localRho(0);

    if(updDens) {

	localRho = localDensity(id);

    }

    else {

	localRho = rho.at(id);

    }


    // Lattice constants

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();

    const uint q = mesh.lmodel()->q();


    // Compute first moment
    
    for( uint j = 0 ; j < 3 ; j++ ) {

    	for( uint k = 0 ; k < q ; k++ ) {

    	    lv[j] += vel[k][j] * _pdf[id][k];
		    
    	}
	    
    }


    // Add interaction force and divide by density

    scalar Ft[3];

    F.total(Ft, id);
    
    for( uint j = 0 ; j < 3 ; j++ )
    	lv[j] = ( lv[j]   +   0.5 * Ft[j]   ) / localRho;
	


    // Copy to global array
    for( uint j = 0 ; j < 3 ; j++ ) 
    	v[j] = lv[j];
	

}







/** Collision process */

const void pseudoPotEquation::collision() {}



/** Update macroscopic density */

const void pseudoPotEquation::updateMacroDensity() {

    for( uint i = 0 ; i < mesh.npoints() ; i++ )
    	rho[i] = localDensity(i);


    
    // for( uint i = 0 ; i < mesh.local() ; i++ )
    // 	rho[i] = localDensity(i);

    // rho.sync();

}


const void pseudoPotEquation::updateMacroDensity(const uint& first, const uint& last) {

    for( uint i = first ; i < last ; i++ )
	rho[i] = localDensity(i);

}




/** Update macroscopic velocity */
/** Repeat for efficiency */

const void pseudoPotEquation::updateMacroVelocity() {

    
    // // Update forces first

    // F.update(rho,T);



    // Lattice constants

    const vector< vector<int> > vel = mesh.lmodel()->lvel();

    const scalar q = mesh.lmodel()->q();

    
    // Local velocity
    
    scalar lv[3] = {0,0,0};

    scalar Ft[3];
    

    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {


    	// Compute first moment
    
    	for( uint j = 0 ; j < 3 ; j++ ) {

    	    lv[j] = 0;
	    
    	    for( uint k = 0 ; k < q ; k++ ) {

    		lv[j] += vel[k][j] * _pdf[i][k];
		    
    	    }
	    
    	}


    	// Add interaction force and divide by density

    	F.total(Ft, i);
    
    	for( uint j = 0 ; j < 3 ; j++ )
    	    U[i][j] = ( lv[j]   +   0.5 * Ft[j]   ) / rho.at(i);

	
    }



    // for( uint i = 0 ; i < mesh.local() ; i++ ) {


    // 	// Compute first moment
    
    // 	for( uint j = 0 ; j < 3 ; j++ ) {

    // 	    lv[j] = 0;
	    
    // 	    for( uint k = 0 ; k < q ; k++ ) {

    // 		lv[j] += vel[k][j] * _pdf[i][k];
		    
    // 	    }
	    
    // 	}


    // 	// Add interaction force and divide by density

    // 	F.total(Ft, i);
    
    // 	for( uint j = 0 ; j < 3 ; j++ )
    // 	    U[i][j] = ( lv[j]   +   0.5 * Ft[j]   ) / localDensity(i);

	
    // }

    // U.sync();

}





/** Forced interface for equilibrium distribution */

const void pseudoPotEquation::eqPS( std::vector<scalar>& n, const scalar& rho_, const std::vector<scalar>& U_ ) const {

    
    // Lattice constants
    
    const uint q = mesh.lmodel()->q();

    const scalar cs2 = mesh.lmodel()->cs2();

    const vector<scalar>& omega = mesh.lmodel()->omega();

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();


    scalar umag(0);

    for( uint j = 0 ; j < 3 ; j++ )
	umag += U_[j] * U_[j];
    
    
    for( uint k = 0 ; k < q ; k++ ) {

	// Dot product

	scalar alpha(0);
	
	for( uint j = 0 ; j < 3 ; j++ )
	    alpha += vel[k][j] * U_[j];

	    
	n[k] = rho_ * omega[k] * ( 1 + alpha/cs2   +   0.5 * alpha * alpha / (cs2*cs2)  -  0.5 * umag / cs2 );

    }

}



/** Update potential as scalar field */

const void pseudoPotEquation::updatePotential( scalarField& phi ) {

    const scalar cs2 = mesh.lmodel()->cs2();
    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

	phi[i] = F.potential( rho.at(i), T.at(i), cs2 );

    }

}



/** Compute and set pressure field */

const void pseudoPotEquation::pressure( const scalarField& phi, scalarField& p ) {

    // Base case. Ideal LBM

    const scalar cs2 = mesh.lmodel()->cs2();
    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

	p[i] = rho.at(i) * cs2;

    }    

}





/** Update contact nodes */

const void pseudoPotEquation::locateContactNodes() {

    
    // Temporary contact nodes
    
    vector<uint> cn;
    
    
    // Move over all boundary nodes and detect density change

    const map< string, vector<uint> >& boundary = mesh.boundaries();

    for( const auto& bd : boundary) {

    	if( bd.second.size() > 0 ) {

    	    for( uint j = 0 ; j < bd.second.size()-1 ; j++ ) {

    		scalar y0( rho.at(bd.second[j]) - _rhoAvgInt ),
    		    y1( rho.at(bd.second[j+1]) - _rhoAvgInt );

    		if( (y1 * y0)   <= 0) {

    		    // cn.push_back( bd.second[j] );
    		    // cn.push_back( bd.second[j+1] );
		    
		    
    		    scalar xc = -y0 / (y1-y0);

    		    uint xint(0);

    		    xc <= 0.5 ?  xint = j : xint = j+1;

		    
    		    if(xint > 0)
    		    	cn.push_back( bd.second[xint-1] );

    		    if(xint - 1 > 0)
    		    	cn.push_back( bd.second[xint-2] );

    		    // if(xint - 2 > 0)
    		    // 	cn.push_back( bd.second[xint-3] );
		    

    		    if(xint + 1 < bd.second.size() )
    		    	cn.push_back( bd.second[xint+1] );

    		    if(xint + 2 < bd.second.size() )
    		    	cn.push_back( bd.second[xint+2] );

    		    // if(xint + 3 < bd.second.size() )
    		    // 	cn.push_back( bd.second[xint+3] );


		    
    		    cn.push_back( bd.second[xint] );
		    

    		}

    	    }

    	}

    }


    _contactNodes = cn;




}
