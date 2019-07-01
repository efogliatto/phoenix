#include "gammaHat.H"

const scalar gammaHat(const latticeMesh& mesh,
		      const scalarField& rho,
		      const scalarField& T,
		      const vectorField& U,
		      const scalar a1,
		      const scalar a2,
		      const scalar tau,
		      const EOS* eos,
		      const scalar Cv,
		      const uint id) {

    
    // Constants

    scalar gradT[3]   = {0,0,0};

    scalar gradRho[3] = {0,0,0};



    // Lattice type

    latticeModel::latticeType ltype = mesh.lmodel()->type();
    


    // Cached scalar values

    const scalar _rho = rho.at(id);

    const scalar _T = T.at(id);
	
	
    
    // Thermal difusivity
    
    scalar chi(0);

    switch( ltype ) {

    case latticeModel::latticeType::D2Q9:

    	chi = (1/tau - 0.5) * (4.0 + 3.0 * a1  + 2.0 * a2) / 6.0;

    	break;


    case latticeModel::latticeType::D3Q15:

    	chi = (1/tau - 0.5) * (6.0 + 11.0 * a1  + a2) / 9.0;

    	break;

    }



    // Scalar fields gradients

    T.grad(gradT, id);

    rho.grad(gradRho, id);

	
    scalar first(0);

    for(uint j = 0 ; j < 3 ; j++)
	first += gradT[j] * gradRho[j];

    first = first * chi / _rho;

	
	
    // Velocity divergence term

    scalar second = U.div(id) * _T * ( 1.0   -   eos->dp_dT(_rho, _T) / (_rho * Cv) );


	
    return first + second;
       
    

}
