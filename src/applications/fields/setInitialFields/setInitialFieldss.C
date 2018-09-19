#include <iostream>

#include <PPEquation.H>


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
    	cout << "     |   - | -   |  " << endl;
    	cout << "     o<----o---->o   Initial fields distribution" << endl;
    	cout << "     |   - | -   |  " << endl;
    	cout << "     | -   |   - |  " << endl;
    	cout << "     o-----o-----o  " << endl << endl;
    }




    // Lattice mesh creation
    
    latticeMesh mesh(pid);


    // Simulation handler

    timeOptions Time(pid);


    // Option dictionary

    dictionary dict("start/initialFieldss");

    
    // Set scalarFields

    vector<string> scfields = dict.bracedEntriesNames("scalarFields");

    for(auto fname : scfields) {


	// Field construction without reading

	scalarField field( mesh, Time, fname, IO::NO_READ, IO::MUST_WRITE );


	// Write field

	field.write();

    }


    // // Macroscopic density

    // scalarField rho( mesh, Time, "rho", IO::MUST_READ, IO::MUST_WRITE );


    // // Macroscopic temperature

    // scalarField T( mesh, Time, "T", IO::MUST_READ, IO::NO_WRITE );


    // // Macroscopic velocity

    // vectorField U( mesh, Time, "U", IO::MUST_READ, IO::MUST_WRITE );
    

    // // PDF field. Navier - Stokes equation

    // pdfField f( mesh, Time, "f", IO::MUST_READ, IO::MUST_WRITE );



   
    
    
    // // Li MRT equation

    // PPEquation NSEq;

    // pseudoPotEquation* NS = NSEq.create("Navier-Stokes", mesh, Time, f, rho, U, T);



    
    // // Advance in time. Collide, stream, update and write
    
    // while( Time.update() ) {

	
	
    // 	// Solve Navier-Stokes equation

    // 	NS->collision();

    // 	NS->streaming();

    // 	NS->updateBoundaries();

    // 	f.sync();


	
    // 	// Update macroscopic fields

    // 	NS->updateMacroDensity();

    // 	NS->updateMacroVelocity();

	

    // 	// Write fields
	
    // 	if( Time.write() ) {


    // 	    rho.write();

    // 	    U.write();

    // 	    f.write();

	    
    // 	    if(pid == 0) {
		
    // 		cout << "Time = " << Time.currentTime() << endl;
		
    // 		cout << "Elapsed time = " << std::fixed << std::setprecision(2) << Time.elapsed() << " seconds" << endl << endl;
		
    // 	    }
	    

    // 	}

    // }






    // // Print info
    // if(pid == 0)	
    // 	cout << endl << "  Finished in " << Time.elapsed() << " seconds " << endl << endl;
	
    

    MPI_Finalize();

}
