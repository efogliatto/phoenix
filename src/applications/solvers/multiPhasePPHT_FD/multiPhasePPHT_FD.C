#include <iostream>

#include <iomanip>

#include <pseudoPotEqHandler.H>

#include <energyEqHandler.H>

#include <armadillo>


using namespace std;

using namespace arma;



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



    // Simulation constants

    dictionary dict("properties/latticeProperties");

    const scalar chi(0.06);

    const uint Nx( dict.lookUp<int>("Nx") );

    const uint Ny( dict.lookUp<int>("Ny") );

    const scalar Tsat(0.094);

    const scalar Tw(0.1367);

    const int spmin(74);

    const int spmax(76);

    const scalar b_eos(2/21);

    const scalar Cv(5);


    
    

    // Energy equation matrix and vector

    const uint armaSize( (Nx * Ny) - 2*Nx );   
    
    arma::vec TVec(armaSize);

    arma::vec newTVec(armaSize);

    arma::vec BVec(armaSize);

    arma::vec K1(armaSize);

    arma::vec K2(armaSize);

    arma::vec K3(armaSize);

    arma::vec K4(armaSize);


       

   

    // Get matrix structure using all neighbours

    vector<int> locX;

    vector<int> locY;

    for( uint i = 0 ; i < armaSize ; i++ ) {

	    
        // LB element id

        int lbid = i + Nx;

        const uint q(9);
	   	   

        for( uint k = 1 ; k < q ; k++ ) {

    	    int nbid = nb[lbid][reverse[k]];

    	    if(  ( (nbid - (int)Nx) > 0 )  &&  ( nbid < (int)(armaSize + Nx) )  ){

		locX.push_back(i);

		locY.push_back(nbid - Nx);

    	    }

        }


    }

    
    umat locMat(2, locX.size());

    vec values(locX.size(), fill::ones);
    
    for(uint i = 0 ; i < locX.size() ; i++) {
	
    	locMat(0,i) = locX[i];

    	locMat(1,i) = locY[i];	

    }


    
    arma::sp_mat TMat(locMat, values);

    arma::sp_mat TMat_1(locMat, values);
    
    arma::sp_mat TMat_2(locMat, values);

    arma::sp_mat TMat_3(locMat, values);

    arma::sp_mat TMat_4(locMat, values);
    
       

    

    // Initialize temperature vector
    
    for( uint i = 0 ; i < armaSize ; i++ )
    	TVec[i] = Tsat;

    
    
    
    // Advance in time. Collide, stream, update and write
    
    while( Time.update() ) {
              
       

    	// Need to initialize coefficients to zero

    	for( vec::iterator it = BVec.begin() ; it != BVec.end() ; ++it )	    	
    	    *it = 0;

	

        // Update TMat	

    	for( uint i = 0 ; i < armaSize ; i++ ) {

	    
    	    // LB element id

    	    int lbid = i + Nx;

    	    const uint q(9);

	    
	    

    	    // Convective term: -U \dot \nabla T

    	    for( uint k = 1 ; k < q ; k++ ) {

    	    	int nbid = nb[lbid][reverse[k]];

    	    	if(  ( (nbid - (int)Nx) > 0 )  &&  ( nbid < (int)(armaSize + Nx) )  ){
    	    	    
    	    	    TMat_1(i, nbid - Nx) = -omega[k] * vel[k][0] * U.at(lbid)[0] / cs2
			                   -omega[k] * vel[k][1] * U.at(lbid)[1] / cs2
    	    		                   -omega[k] * vel[k][2] * U.at(lbid)[2] / cs2;

    	    	}

    	    	else {

    	    	    scalar Tb(Tsat);

		    if( mesh.latticePoint(nbid)[1] == 0  ) {
		    
			if(   ( mesh.latticePoint(nbid)[0] >= spmin )  &&  ( mesh.latticePoint(nbid)[0] <= spmax )  ) {
			    
			    Tb = Tw;

			}

		    }

    	    	    for( uint j = 0 ; j < 3 ; j++ )		    
    	    	    	BVec(i) -= Tb * omega[k] * vel[k][j] * U.at(lbid)[j] / cs2;

    	    	}

    	    }





    	    // Diffusive term: chi \nabla^2 T

	    TMat_2.at(i, i) = 0;

    	    for( uint k = 1 ; k < q ; k++ ) {

    	    	int nbid = nb[lbid][reverse[k]];

    	    	if(  ( (nbid - (int)Nx) > 0 )  &&  ( nbid < (int)(armaSize + Nx) )  ){

    	    	    TMat_2.at(i, nbid - Nx) = 2 * omega[k] * chi / cs2;

    	    	    TMat_2.at(i, i) -= 2 * omega[k] * chi / cs2;		    

    	    	}

    	    	else {

    	    	    scalar Tb(Tsat);

		    if( mesh.latticePoint(nbid)[1] == 0 ) {
		    
			if(   ( mesh.latticePoint(nbid)[0] >= spmin )  &&  ( mesh.latticePoint(nbid)[0] <= spmax ) ) {
			    
			    Tb = Tw;

			}

		    }

    	    	    BVec(i) += 2 * Tb * omega[k] * chi / cs2;

		    TMat_2.at(i, i) -= 2 * omega[k] * chi / cs2;

    	    	}

    	    }





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

    	    // 	if(  ( (nbid - (int)Nx) > 0 )  &&  ( nbid < (int)(armaSize + Nx) )  ){
    	    	   
    	    // 	    TMat_3.at(i, nbid - Nx) = chi * omega[k] * vel[k][0] * gradRho[0] / (cs2 * rho.at(lbid))
    	    // 		                    + chi * omega[k] * vel[k][1] * gradRho[1] / (cs2 * rho.at(lbid))
    	    // 		                    + chi * omega[k] * vel[k][2] * gradRho[2] / (cs2 * rho.at(lbid));

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
	    

    	    // TMat_4(i,i) = -dpdT * divU / ( rho.at(lbid) * Cv );
	   	    
	    

    	}



	

    	// Construct RK vectors and advance in time

    	TMat = TMat_1 + TMat_2;// + TMat_3 + TMat_4;
	
    	K1 = (TMat * TVec) + BVec;

    	K2 = TMat * ( TVec + 0.5*K1 ) + BVec; 

    	K3 = TMat * ( TVec + 0.5*K2 ) + BVec; 

    	K4 = TMat * ( TVec + K3 )     + BVec;

    	TVec = TVec + (1.0/6.0)*K1 + (1.0/3.0)*(K2+K3) + (1.0/6.0)*K4;

    	// TVec += newTVec;


	
	
    	// Update T field

    	for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

    	    if( mesh.latticePoint(i)[1] == 0 ) {

    		if(   ( mesh.latticePoint(i)[0] >= spmin )  &&  ( mesh.latticePoint(i)[0] <= spmax )  ) {

    		    T[i] = Tw;

    		}

    		else {

    		    T[i] = Tsat;

    		}
			

    	    }

    	    else {

    		if( mesh.latticePoint(i)[1] == (int)(Ny - 1) ) {

    		    T[i] = Tsat;

    		}

    		else {

    		    T[i] = TVec.at(i-Nx);

    		}

    	    }

	    

    	}
    





    	// Solve Navier-Stokes equation

    	NS.collision();

    	NS.streaming();

    	NS.updateBoundaries();

    	f.sync();

    	NS.updateMacroDensity();

    	NS.updateMacroVelocity();


	

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
