#include <iostream>

#include <iomanip>

#include <pseudoPotEqHandler.H>

#include <energyEqHandler.H>

#include <armadillo>


using namespace std;



int main( int argc, char **argv ) {

    
    // Initialize mpi

    int pid, world;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD,&pid);

    MPI_Comm_size(MPI_COMM_WORLD,&world);



    if(pid == 0) {
    	cout << "                    " << endl;
    	cout << "     o-----o-----o  " << endl;
    	cout << "     | -   |   - |  " << endl;
    	cout << "     |   - | -   |  Two Phases - Lattice-Boltzmann solver with heat transfer" << endl;
    	cout << "     o<----o---->o  " << endl;
    	cout << "     |   - | -   |                Pseudopotential model" << endl;
    	cout << "     | -   |   - |  " << endl;
    	cout << "     o-----o-----o  " << endl << endl;
    }




    // Lattice mesh creation
    
    latticeMesh mesh(pid);


    // Simulation handler

    timeOptions Time(pid);


    // Macroscopic density

    scalarField rho( mesh, Time, "rho", IO::MUST_READ, IO::MUST_WRITE );


    // Macroscopic temperature

    scalarField T( mesh, Time, "T", IO::MUST_READ, IO::MUST_WRITE );


    // Macroscopic velocity

    vectorField U( mesh, Time, "U", IO::MUST_READ, IO::MUST_WRITE );
    

    // PDF field. Navier - Stokes equation

    pdfField f( mesh, Time, "f", IO::MUST_READ, IO::MUST_WRITE );


   
    
    
    // Navier-Stokes MRT equation

    pseudoPotEqHandler NS("Navier-Stokes", mesh, Time, f, rho, U, T);




    // Mesh neighbours

    const vector< vector<int> >& nb = mesh.nbArray();

    const vector<uint>& reverse = mesh.lmodel()->reverse();

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();

    const scalar cs2 = mesh.lmodel()->cs2();

    const vector<scalar>& omega = mesh.lmodel()->omega();    


    

    // Energy equation matrix and vector

    const uint armaSize(45149);
        
    arma::sp_mat TMat(armaSize, armaSize);   

    arma::vec TVec(armaSize);

    arma::vec BVec(armaSize);

    arma::vec K1(armaSize);

    arma::vec K2(armaSize);

    arma::vec K3(armaSize);

    arma::vec K4(armaSize);


       

   
    // Simulation constants

    const scalar chi(0.06);

    const uint Nx(151);

    const scalar Tsat(0.094);

    const scalar Tw(0.1367);

    const int spmin(74);

    const int spmax(76);

    const scalar b_eos(2/21);

    const scalar Cv(5);    

    

    // Initialize temperature vector
    
    for( uint i = 0 ; i < armaSize ; i++ )
	TVec[i] = Tsat;
    
    
    
    // Advance in time. Collide, stream, update and write
    
    while( Time.update() ) {

	
       
	
	// Solve Navier-Stokes equation

	NS.collision();

	NS.streaming();

	NS.updateBoundaries();

	f.sync();

	NS.updateMacroDensity();

	NS.updateMacroVelocity();
	


	// Update TMat

	// Need to initialize coefficients to zero

	for( uint i = 0 ; i < armaSize ; i++ ) {

	    
	    // LB element id

	    int lbid = i + 151;

	    const uint q(9);

	    
	    

	    // Convective term: -U \dot \nabla T

	    for( uint k = 1 ; k < q ; k++ ) {

		int nbid = nb[lbid][reverse[k]];

		if(  ( (nbid - Nx) > 0 )  &&  ( nbid < (int)armaSize + (int)Nx )  ){

		    // for( uint j = 0 ; j < 3 ; j++ )
		    	// TMat(i, nbid - Nx) -= omega[k] * vel[k][j] * U.at(lbid)[j] / cs2;

		}

		else {

		    scalar Tb(Tsat);
		    
		    if(   ( mesh.latticePoint(lbid)[0] >= spmin )  &&  ( mesh.latticePoint(lbid)[0] <= spmax )  )
		    	Tb = Tw;

		    for( uint j = 0 ; j < 3 ; j++ )		    
		    	BVec(i) -= Tb * omega[k] * vel[k][j] * U.at(lbid)[j] / cs2;

		}

	    }





	    // // Diffusive term: chi \nabla^2 T

	    // for( uint k = 1 ; k < q ; k++ ) {

	    // 	int nbid = nb[lbid][reverse[k]];

	    // 	if(  ( (nbid - Nx) > 0 )  &&  ( (uint)nbid < armaSize + Nx )  ){

	    // 	    TMat(i, nbid - Nx) += 2 * omega[k] * chi / cs2;

	    // 	    TMat(i, i) -= 2 * omega[k] * chi / cs2;		    

	    // 	}

	    // 	else {

	    // 	    scalar Tb(Tsat);

	    // 	    if(   ( mesh.latticePoint(lbid)[0] >= spmin )  &&  ( mesh.latticePoint(lbid)[0] <= spmax )  )
	    // 	    	Tb = Tw;

	    // 	    BVec(i) += 2 * Tb * omega[k] * chi / cs2;

	    // 	}

	    // }





	    // // Diffusive term: chi (\nabla rho) \cdot (\nabla T) / \rho


	    // // Density gradient
	    
	    // scalar gradRho[3] = {0,0,0};

	    // for( uint k = 1 ; k < q ; k++ ) {

	    // 	int nbid = nb[lbid][reverse[k]];

	    // 	for( uint j = 0 ; j < 3 ; j++ )
	    // 	    gradRho[j] += omega[k] * rho.at(nbid) * vel[k][j] / cs2; 
		

	    // }
	    
	   
	    // for( uint k = 1 ; k < q ; k++ ) {

	    // 	int nbid = nb[lbid][reverse[k]];

	    // 	if(  ( (nbid - Nx) > 0 )  &&  ( (uint)nbid < armaSize + Nx )  ){

	    // 	    for( uint j = 0 ; j < 3 ; j++ )
	    // 	    	TMat(i, nbid - Nx) += chi * omega[k] * vel[k][j] * gradRho[j] / (cs2 * rho.at(lbid));

	    // 	}

	    // 	else {

	    // 	    scalar Tb(Tsat);

	    // 	    if(   ( mesh.latticePoint(lbid)[0] >= spmin )  &&  ( mesh.latticePoint(lbid)[0] <= spmax )  )
	    // 	    	Tb = Tw;

	    // 	    for( uint j = 0 ; j < 3 ; j++ )		    
	    // 	    	BVec(i) += chi * Tb * omega[k] * vel[k][j] * gradRho[j] /  (cs2 * rho.at(lbid));

	    // 	}

	    // }






	    // // Extra term

	    // scalar divU = 0.5 * U.at( nb[lbid][3] )[0]
	    // 	        - 0.5 * U.at( nb[lbid][1] )[0]
	    // 	        + 0.5 * U.at( nb[lbid][4] )[1]
	    // 	        - 0.5 * U.at( nb[lbid][2] )[1];

	    // scalar dpdT( rho.at(lbid) / (1 - rho.at(lbid) * b_eos) );
	    

	    // TMat(i,i) -= dpdT * divU / ( rho.at(lbid) * Cv );
	   	    

	    

	}




	


	// // Construct RK vectors and advance in time

	// K1 = TMat * TVec + BVec;

	// K2 = TMat * ( TVec + 0.5*K1 ) + BVec;

	// K3 = TMat * ( TVec + 0.5*K2 ) + BVec;

	// K4 = TMat * ( TVec + K3 )     + BVec;

	// TVec = (1/6)*K1 + (1/3)*(K2+K3) + (1/6)*K4;

	
	
	


	

    	// Write fields
	
    	if( Time.write() ) {


    	    rho.write();

    	    U.write();

	    T.write();

    	    f.write();


	    
    	    if(pid == 0) {
		
    		cout << "Time = " << Time.currentTime() << endl;
		
    		cout << "Elapsed time = " << std::fixed << std::setprecision(2) << Time.elapsed() << " seconds" << endl << endl;
		
    	    }
	    

    	}

    }






    // Print info
    if(pid == 0)	
    	cout << endl << "  Finished in " << Time.elapsed() << " seconds " << endl << endl;
	
    

    MPI_Finalize();

}
