#include <myMRTEq.H>

#include <sparseScalarMatrix.H>

using namespace std;


/** Default constructor */

myMRTEq::myMRTEq( const std::string& name,
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
    


    // Set initial values to equilibrium distribution

    myMRTEq::setEquilibrium();

}



/** Default destructor */

myMRTEq::~myMRTEq() {}



/** Equilibrium in moment space */

const void myMRTEq::eqMS( vector<scalar>& m, const uint& id ) const {

    
    const uint q = mesh.lmodel()->q();

    const scalar _T = T.at(id);

    const scalar _U[3] = { U.at(id,0), U.at(id,1), U.at(id,2) };
    

    switch(q) {

    case 9:

	m[0] = _T;
	
    	m[1] = _a1 * _T;
	
    	m[2] = _a2 * _T;
	
    	m[3] = _T * _U[0];
	
    	m[4] = _T * (-_U[0]);
	
    	m[5] = _T * _U[1];
	
    	m[6] = _T * (-_U[1]);
	
    	m[7] = 0;
	
    	m[8] = 0;


    	break;


    case 15:

	m[0]  = _T;
	
	m[1]  = _a1 * _T;

	m[2]  = _a2 * _T;

	m[3]  =   _T * _U[0];

	m[4]  = -(_T * _U[0]);

	m[5]  =   _T * _U[1];

	m[6]  = -(_T * _U[1]);

	m[7]  =   _T * _U[2];

	m[8]  = -(_T * _U[2]);

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

    }

}





/** Equilibrium in population space */

const void myMRTEq::eqPS( vector<scalar>& n, const uint& id ) const {


    // First compute in moment space

    myMRTEq::eqMS(n,id);


    
    // Back to population space

    const uint q = mesh.lmodel()->q();
    
    const scalarMatrix& invM = mesh.lmodel()->MRTInvMatrix();

    vector<scalar> res(q);

    invM.matDotVec(n,res);


    std::copy(res.begin(), res.end(), n.begin());


}


/** Set pdf to equilibrium values */

const void myMRTEq::setEquilibrium() {

    const scalar q = mesh.lmodel()->q();
    
    vector<scalar> n( q );
    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

	eqPS(n,i);

	for(uint k = 0 ; k < q ; k++)
	    _pdf.set(i,k,n[k]);

    }

}



/** Collision process */

const void myMRTEq::collision() {


    
    // Lattice constants

    scalarMatrix M = mesh.lmodel()->MRTMatrix();

    scalarMatrix invM = mesh.lmodel()->MRTInvMatrix();    

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

	cout << " [ERROR]  Undefined grid model for myMRT" << endl << endl;

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

const void myMRTEq::updateMacroTemperature() {

    // Update heat sources first

    _hs->update(rho,T,U);

    
    for( uint i = 0 ; i < mesh.npoints() ; i++ )
	T[i] = myMRTEq::localTemperature(i);

}
