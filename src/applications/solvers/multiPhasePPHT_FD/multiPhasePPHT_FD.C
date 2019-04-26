#include <iostream>

#include <iomanip>

#include <pseudoPotEqHandler.H>

#include <energyEqHandler.H>

#include <TEquation.H>


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


    // Finite-difference T equation

    TEquation Teq( mesh, Time, T);    


    // Runge-Kutta fields

    scalarField K0( mesh, Time, "K0", IO::NO_READ, IO::NO_WRITE );

    scalarField K1( mesh, Time, "K1", IO::NO_READ, IO::NO_WRITE );

    scalarField K2( mesh, Time, "K2", IO::NO_READ, IO::NO_WRITE );

    scalarField K3( mesh, Time, "K3", IO::NO_READ, IO::NO_WRITE );

    scalarField K4( mesh, Time, "K4", IO::NO_READ, IO::NO_WRITE );        
    



    // PDF field. Energy equation

    pdfField g( mesh, Time, "g", IO::MUST_READ, IO::MUST_WRITE );


    // Macroscopic temperature

    scalarField Tlb( mesh, Time, "Tlb", IO::NO_READ, IO::MUST_WRITE );

    for(uint i = 0 ; i < mesh.npoints() ; i++)
	Tlb[i] = T.at(i);


    
    // Energy MRT equation

    energyEqHandler energy("Energy", mesh, Time, g, rho, U, Tlb);
    
    
    
    
    // Advance in time. Collide, stream, update and write
    
    while( Time.update() ) {
              
       
	// // Energy equation. RK-4 scheme

	
	// // K1
	
	// Teq.rval(K1, rho, U, T);


	// // K2

	// for( uint id = 0 ; id < mesh.local() ; id++ )
	//     K0[id] = T.at(id) + 0.5*K1.at(id);

	// K0.sync();

	// Teq.rval(K2, rho, U, K0);


	// // K3

	// for( uint id = 0 ; id < mesh.local() ; id++ )
	//     K0[id] = T.at(id) + 0.5*K2.at(id);

	// K0.sync();

	// Teq.rval(K3, rho, U, K0);


	// // K4

	// for( uint id = 0 ; id < mesh.local() ; id++ )
	//     K0[id] = T.at(id) + K3.at(id);

	// K0.sync();

	// Teq.rval(K4, rho, U, K0);



	// // Update T

	// for( uint id = 0 ; id < mesh.local() ; id++ )
	//     T[id] = T.at(id) + (1.0/6.0)*K1.at(id) + (1.0/3.0)*K2.at(id) + (1.0/3.0)*K3.at(id) + (1.0/6.0)*K4.at(id);

	// T.sync();




	


	// Energy equation. RK-2 scheme

	
	// K1
	
	Teq.rval(K1, rho, U, T);


	// K2

	for( uint id = 0 ; id < mesh.local() ; id++ )
	    K0[id] = T.at(id) + 0.5*K1.at(id);

	K0.sync();

	Teq.rval(K2, rho, U, K0);


	// Update T

	for( uint id = 0 ; id < mesh.local() ; id++ )
	    T[id] = T.at(id) + (1.0/2.0)*K1.at(id) + (1.0/2.0)*K2.at(id);

	T.sync();



	

	// // Energy equation. Explicit Euler scheme

	
	// // K1
	
	// Teq.rval(K1, rho, U, T);


	// // Update T

	// for( uint id = 0 ; id < mesh.local() ; id++ )
	//     T[id] = T.at(id) + K1.at(id);

	// T.sync();

	


	// Solve Energy equation

	energy.collision();

	energy.streaming();

	energy.updateBoundaries();

	g.sync();

	energy.updateMacroTemperature();

	for(uint i = 0 ; i < mesh.local() ; i++) {

	    if( mesh.isOnBoundary(i) ) {

		Tlb[i] = T.at(i);

	    }

	}

	T.sync();

	
	

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

	    Tlb.write();

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
