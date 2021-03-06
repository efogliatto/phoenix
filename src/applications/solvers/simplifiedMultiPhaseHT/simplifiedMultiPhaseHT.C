#include <iostream>

#include <iomanip>

#include <pseudoPotEqHandler.H>

#include <energyEqHandler.H>

#include <simplifiedTEq.H>



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
    	cout << "     |   - | -   |  Two Phases - Lattice Boltzmann solver with heat transfer" << endl;
    	cout << "     o<----o---->o  " << endl;
    	cout << "     |   - | -   |                Pseudopotential model" << endl;
    	cout << "     | -   |   - |            Predictor - corrector scheme" << endl;
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


    // PDF field. Energy equation

    pdfField g( mesh, Time, "g", IO::MUST_READ, IO::MUST_WRITE );    
    





    // Finite-difference T equation

    scalarField Ts( mesh, Time, "Ts", IO::MUST_READ, IO::MUST_WRITE );

    scalarField Tstar( mesh, Time, "Tstar", IO::NO_READ, IO::NO_WRITE );    
    
    simplifiedTEq Teq( mesh, Time, Ts );    


    // Navier-Stokes MRT equation

    // pseudoPotEqHandler NS("Navier-Stokes", mesh, Time, f, rho, U, T);
    pseudoPotEqHandler NS("Navier-Stokes", mesh, Time, f, rho, U, Ts);    


    // // Energy MRT equation

    // energyEqHandler energy("Energy", mesh, Time, g, rho, U, T);    
 

    

    
    
    // Advance in time. Collide, stream, update and write
    
    while( Time.update() ) {
              


	// Energy equation

	
	// Predictor step

	Teq.predictor(Tstar, rho, U, Ts);

	    
	// Corrector step

	Teq.corrector(Ts, rho, U, Tstar);
	

	



	// // Solve Energy equation

	// energy.collision();

	// energy.streaming();

	// energy.updateBoundaries();

	// g.sync();

	// energy.updateMacroTemperature();

	

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

	    Ts.write();

    	    f.write();

    	    g.write();	    


	    
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
