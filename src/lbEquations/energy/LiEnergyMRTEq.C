#include <LiEnergyMRTEq.H>

#include <algebra.H>

using namespace std;


/** Default constructor */

LiEnergyMRTEq::LiEnergyMRTEq( const std::string& name,
			      const latticeMesh& mesh_,
			      timeOptions& Time_,
			      pdfField& pdf_,
			      const scalarField& rho_,
			      const vectorField& U_,
			      scalarField& T_)

    : energyEquation(name, mesh_, Time_, pdf_, rho_, U_, T_) {


    // Read constants

    dictionary dict("properties/macroProperties");

    _Cv = dict.lookUp<scalar>( name + "/HeatSource/Constants/Cv");

    _kappa = dict.lookUp<scalar>( name + "/HeatSource/Constants/kappa");
    


    // Set equilibrium pdf values only if time is start time

    if( Time.currentTime() == Time.timeList()[0] )  {

	LiEnergyMRTEq::setEquilibrium();

    }

}



/** Default destructor */

LiEnergyMRTEq::~LiEnergyMRTEq() {}



/** Equilibrium in moment space */

const void LiEnergyMRTEq::eqMS( vector<scalar>& m, const uint& id ) const {
  
    LiEnergyMRTEq::eqMS( m, T.at(id), U.at(id) );

}



/** Equilibrium in moment space */

const void LiEnergyMRTEq::eqMS( vector<scalar>& m, const scalar& T_, const vector<scalar>& U_ ) const {
   

    const scalar _U[3] = { U_[0], U_[1], U_[2] };
    

    switch( mesh.lmodel()->type()) {

    case latticeModel::latticeType::D2Q9:

	m[0] = T_;
	
    	m[1] = -2 * T_;
	
    	m[2] = 2 * T_;
	
    	m[3] = T_ * _U[0];
	
    	m[4] = T_ * (-_U[0]);
	
    	m[5] = T_ * _U[1];
	
    	m[6] = T_ * (-_U[1]);
	
    	m[7] = 0;
	
    	m[8] = 0;


    	break;



    default:

    	cout << " [ERROR]  Equilibrium model not implemented" << endl;

    	exit(1);

	break;

    }

}






/** Equilibrium in population space */

const void LiEnergyMRTEq::eqPS( vector<scalar>& n, const uint& id ) const {

    LiEnergyMRTEq::eqPS( n, T.at(id), U.at(id), 0 );
    
}



/** Equilibrium in population space */

const void LiEnergyMRTEq::eqPS( vector<scalar>& n, const scalar& T_, const std::vector<scalar>& U_, const scalar& hs ) const {



    // First compute in moment space
   
    vector<scalar> n_eq( mesh.lmodel()->q() );

    LiEnergyMRTEq::eqMS(n_eq,T_,U_);


    
    // Back to population space
    
    const scalarMatrix& invM = mesh.lmodel()->MRTInvMatrix();

    invM.matDotVec(n_eq,n);

}



/** Set pdf to equilibrium values */

const void LiEnergyMRTEq::setEquilibrium() {

    const scalar q = mesh.lmodel()->q();
    
    vector<scalar> n( q );
    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

	eqPS(n,i);

	for(uint k = 0 ; k < q ; k++) 
	    _pdf.set(i,k,n[k]);	    

    }

}



/** Collision process */

const void LiEnergyMRTEq::collision() {


    
    // Lattice constants

    const scalarMatrix& M = mesh.lmodel()->MRTMatrix();

    const scalarMatrix& invM = mesh.lmodel()->MRTInvMatrix();    

    const uint nodes = mesh.npoints();

    const uint q = mesh.lmodel()->q();


    
    // Partial distributions
    
    vector<scalar> n(q);     // m:  momentum space
    
    vector<scalar> n_eq(q);  // meq: equilibrium in momentum space

    
    

    // Move over all points

    for( uint id = 0 ; id < nodes ; id++ ) {

	
    	// Compute equilibrium in moment space

    	eqMS(n_eq,id);


    	// Distribution in moment space

    	M.matDotVec( _pdf[id], n );

	scalar n4 = n[4];

	scalar n6 = n[6];


    	// Collision in momentum space
	
    	for( uint k = 0 ; k < q ; k++ )
    	    n[k] = n[k]   -   _Tau[k] * (n[k] - n_eq[k]);

	    
    	n[0] = n[0] + _hs->source(id);



	// Correction terms

	n[3] = n[3]   +   ( 1 - 0.5*_Tau[3] ) * _Tau[4] * (n4 - n_eq[4]);

	n[5] = n[5]   +   ( 1 - 0.5*_Tau[5] ) * _Tau[6] * (n6 - n_eq[6]);
	


    	// Back to population space

    	invM.matDotVec(n, _pdf[id]);
	

    }



}







/** Update macroscopic temperature */

const void LiEnergyMRTEq::updateMacroTemperature() {

    
    // Update heat sources first

    _hs->update(rho,T,U);

    
    for( uint i = 0 ; i < mesh.npoints() ; i++ )
	T[i] = LiEnergyMRTEq::localTemperature(i);


    // _hs->update(rho,T,U);    
}




/** Thermal conductivity at node */

const scalar LiEnergyMRTEq::thermalCond( const uint& id ) const {

    cout << "Thermal cond not fully verified" << endl;

    exit(1);
    
    return _kappa;

}
