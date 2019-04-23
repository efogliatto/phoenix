#include <GongSRTEq.H>

#include <algebra.H>

using namespace std;


/** Default constructor */

GongSRTEq::GongSRTEq( const std::string& name,
		  const latticeMesh& mesh_,
		  timeOptions& Time_,
		  pdfField& pdf_,
		  const scalarField& rho_,
		  const vectorField& U_,
		  scalarField& T_)

    : energyEquation(name, mesh_, Time_, pdf_, rho_, U_, T_) {


    // Read constants

    dictionary dict("properties/macroProperties");

    _a1 = dict.lookUp<scalar>( name + "/HeatSource/Constants/alpha_1");

    _a2 = dict.lookUp<scalar>( name + "/HeatSource/Constants/alpha_2");

    _Cv = dict.lookUp<scalar>( name + "/HeatSource/Constants/Cv");
    


    // Set equilibrium pdf values only if time is start time

    if( Time.currentTime() == Time.timeList()[0] )  {

	GongSRTEq::setEquilibrium();

    }

}



/** Default destructor */

GongSRTEq::~GongSRTEq() {}



/** Equilibrium in moment space */

const void GongSRTEq::eqMS( vector<scalar>& m, const uint& id ) const {
  
    GongSRTEq::eqMS( m, T.at(id), U.at(id) );

}



/** Equilibrium in moment space */

const void GongSRTEq::eqMS( vector<scalar>& m, const scalar& T_, const vector<scalar>& U_ ) const {

    
    const uint q = mesh.lmodel()->q();

    const scalar _U[3] = { U_[0], U_[1], U_[2] };
    

    switch(q) {

    case 9:

	m[0] = T_;
	
    	m[1] = _a1 * T_;
	
    	m[2] = _a2 * T_;
	
    	m[3] = T_ * _U[0];
	
    	m[4] = T_ * (-_U[0]);
	
    	m[5] = T_ * _U[1];
	
    	m[6] = T_ * (-_U[1]);
	
    	m[7] = 0;
	
    	m[8] = 0;


    	break;


    case 15:

	m[0]  = T_;
	
	m[1]  = _a1 * T_;

	m[2]  = _a2 * T_;

	m[3]  =   T_ * _U[0];

	m[4]  = -(T_ * _U[0]);

	m[5]  =   T_ * _U[1];

	m[6]  = -(T_ * _U[1]);

	m[7]  =   T_ * _U[2];

	m[8]  = -(T_ * _U[2]);

	m[9]  = 0;

	m[10] = 0;

	m[11] = 0;

	m[12] = 0;

	m[13] = 0;

	m[14] = 0;
	
    	break;



    default:

    	cout << " [ERROR]  Equilibrium model not implemented" << endl;

    	exit(1);

	break;

    }

}






/** Equilibrium in population space */

const void GongSRTEq::eqPS( vector<scalar>& n, const uint& id ) const {

    GongSRTEq::eqPS( n, T.at(id), U.at(id), _hs->source(id) );
    
}



/** Equilibrium in population space */

const void GongSRTEq::eqPS( vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_, const scalar& hs ) const {


    const uint q = mesh.lmodel()->q();
    

    // First compute in moment space
   
    vector<scalar> n_eq( q );

    GongSRTEq::eqMS(n_eq,T_,U_);


    
    // Back to population space
    
    const scalarMatrix& invM = mesh.lmodel()->MRTInvMatrix();

    invM.matDotVec(n_eq,n);



    // Source term in moment space
   
    vector<scalar> GammaHat(q);

    vector<scalar> Gamma(q);

    GammaHat[0] = hs;

    invM.matDotVec(GammaHat,Gamma);
    

    
    // Update Source term

    for( uint k = 0 ; k < q ; k++ )    
	n[k] = n[k] - 0.5 * Gamma[k];

}



/** Set pdf to equilibrium values */

const void GongSRTEq::setEquilibrium() {

    const scalar q = mesh.lmodel()->q();
    
    vector<scalar> n( q );
    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

	eqPS(n,i);

	for(uint k = 0 ; k < q ; k++) 
	    _pdf.set(i,k,n[k]);	    

    }

}



/** Collision process */

const void GongSRTEq::collision() {


    
    // Lattice constants

    const scalarMatrix& M = mesh.lmodel()->MRTMatrix();

    const scalarMatrix& invM = mesh.lmodel()->MRTInvMatrix();    

    const uint nodes = mesh.npoints();

    const uint q = mesh.lmodel()->q();


    
    // Partial distributions
    
    vector<scalar> n(q);     // m:  momentum space
    
    vector<scalar> n_eq(q);  // meq: equilibrium in momentum space

    vector<scalar> GammaHat(q);     // Source in population space

    vector<scalar> aux_1(q);

    vector<scalar> aux_2(q);
    


    // Non diagonal Q

    sparseScalarMatrix Q( _Tau );

    switch(q) {

    case 9:

    	Q.addElement( _Tau[4]  *  ( _Tau[3]/2.0  - 1.0 ), 3, 4);

    	Q.addElement( _Tau[6]  *  ( _Tau[5]/2.0  - 1.0 ), 5, 6);
	
    	break;


    case 15:

    	Q.addElement( _Tau[4]  *  ( _Tau[3]/2.0  - 1.0 ), 3, 4);

    	Q.addElement( _Tau[6]  *  ( _Tau[5]/2.0  - 1.0 ), 5, 6);

    	Q.addElement( _Tau[8]  *  ( _Tau[7]/2.0  - 1.0 ), 7, 8);	
	
    	break;

    default:

	cout << " [ERROR]  Undefined grid model for GongSRT" << endl << endl;

	exit(1);
	
	break;

    }

    

    // Move over all points

    for( uint id = 0 ; id < nodes ; id++ ) {

	
    	// Compute equilibrium in moment space

    	eqMS(n_eq,id);


    	// Distribution in moment space

    	M.matDotVec( _pdf[id], n );



    	// First auxiliary distribution: n - n_eq

    	for( uint k = 0 ; k < q ; k++ )
    	    aux_2[k] = n[k] - n_eq[k];


	
    	// aux_1 = Q * aux_2 = Q * (n - n_eq)

    	Q.matDotVec( aux_2, aux_1 );


	
        // Second auxiliary distribution: (I  -  0.5 * Q) * GammaHat	

    	scalar heat = (1.0 - 0.5 * _Tau[0])  * _hs->source(id);



    	// Collision in momentum space
	
    	for( uint k = 0 ; k < q ; k++ )
    	    n[k] = n[k] - aux_1[k];

	
    	n[0] += heat;



    	// Back to population space

    	invM.matDotVec(n, _pdf[id]);
	

    }



}







/** Update macroscopic temperature */

const void GongSRTEq::updateMacroTemperature() {

    // Update heat sources first

    _hs->update(rho,T,U);

    
    for( uint i = 0 ; i < mesh.npoints() ; i++ )
	T[i] = GongSRTEq::localTemperature(i);

}




/** Thermal conductivity at node */

const scalar GongSRTEq::thermalCond( const uint& id ) const {

    scalar k(0);
    
    
    switch( mesh.lmodel()->q()) {

    case 9:

	k = rho.at(id) * _Cv * (1/_Tau[3] - 0.5) * (4.0 + 3.0 * _a1  + 2.0 * _a2) / 6.0;
	
	break;


    case 15:

	k = rho.at(id) * _Cv * (1/_Tau[3] - 0.5) * (6.0 + 11.0 * _a1  + _a2) / 9.0;
	    	
	break;


    }


    return k;

}
