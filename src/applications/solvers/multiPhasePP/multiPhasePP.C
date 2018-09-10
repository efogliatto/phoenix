#include <iostream>

#include <iomanip>

#include <scalarField.H>

#include <vectorField.H>

#include <pdfField.H>


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

    timeOptions Time;


    // Macroscopic density

    scalarField rho( mesh, Time, "rho" );


    // Macroscopic temperature

    scalarField T( mesh, Time, "T" );


    // Macroscopic velocity

    vectorField U( mesh, Time, "U" );
    

    // PDF field. Navier - Stokes equation

    pdfField f( mesh, Time, "f" );


    // // PDF field. energy equation

    // pdfField g( mesh, Time, "g" );

    



    // Advance in time. Collide, stream, update and write
    
    while( Time.update() ) {



	// Write fields
	
	if( Time.write() ) {


	    rho.write();

	    T.write();

	    U.write();

	    f.write();

	    // g.write();

	    
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
