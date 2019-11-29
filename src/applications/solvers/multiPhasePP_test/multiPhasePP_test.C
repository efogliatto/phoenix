#include <iostream>

#include <iomanip>

#include <pseudoPotEqHandler.H>

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
    	cout << "     |   - | -   |  Two Phases - Lattice-Boltzmann solver" << endl;
    	cout << "     o<----o---->o  " << endl;
    	cout << "     |   - | -   |          Pseudopotential model" << endl;
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




    // Timing
    
    typedef std::chrono::high_resolution_clock clock;

    typedef std::chrono::duration<double, std::ratio<1> > second;

    std::chrono::time_point<clock> beg;

    vector<scalar> elapsed(6);
    std::fill( elapsed.begin(), elapsed.end(), 0 );
    

    
   
    
    
    // Li MRT equation

    pseudoPotEqHandler NS("Navier-Stokes", mesh, Time, f, rho, U, T);



    
    // Advance in time. Collide, stream, update and write
    
    while( Time.update() ) {

	
	
	// Solve Navier-Stokes equation

	beg = clock::now();

	NS.collision();

	elapsed[0] += (scalar)chrono::duration_cast<second> ( clock::now() - beg ).count();


	
	beg = clock::now();

	NS.streaming();

	elapsed[1] += (scalar)chrono::duration_cast<second> ( clock::now() - beg ).count();	

	

	beg = clock::now();

	NS.updateBoundaries();

	elapsed[2] += (scalar)chrono::duration_cast<second> ( clock::now() - beg ).count();
	


	beg = clock::now();

	f.sync();

	elapsed[3] += (scalar)chrono::duration_cast<second> ( clock::now() - beg ).count();

	


	
	// Update macroscopic fields

	beg = clock::now();
	
	NS.updateMacroDensity();

	elapsed[4] += (scalar)chrono::duration_cast<second> ( clock::now() - beg ).count();

	


	beg = clock::now();
	
	NS.updateMacroVelocity();

	elapsed[5] += (scalar)chrono::duration_cast<second> ( clock::now() - beg ).count();	

	

    	// Write fields
	
    	if( Time.write() ) {


    	    rho.write();

    	    U.write();

    	    f.write();

	    T.write();

	    
    	    if(pid == 0) {
		
    		cout << "Time = " << Time.currentTime() << endl;
		
    		cout << "Elapsed time = " << std::fixed << std::setprecision(2) << Time.elapsed() << " seconds" << endl << endl;
		
    	    }
	    

    	}

    }






    // Print info
    if(pid == 0) {
	
    	cout << endl << "  Finished in " << Time.elapsed() << " seconds " << endl << endl;

	cout << "  Collision " << elapsed[0] << " seconds " << endl;
	cout << "  Streaming " << elapsed[1] << " seconds " << endl;
	cout << "  Boundary  " << elapsed[2] << " seconds " << endl;
	cout << "  Sync      " << elapsed[3] << " seconds " << endl;
	cout << "  Density   " << elapsed[4] << " seconds " << endl;
	cout << "  Velocity  " << elapsed[5] << " seconds " << endl << endl;	

    }
	
    

    MPI_Finalize();

}
