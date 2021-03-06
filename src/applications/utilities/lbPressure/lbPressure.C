#include <iostream>

#include <iomanip>

#include <pseudoPotEqHandler.H>

#include <energyEqHandler.H>


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
    	cout << "     o<----o---->o  Pressure calculation" << endl;
    	cout << "     |   - | -   |  " << endl;
    	cout << "     | -   |   - |  " << endl;
    	cout << "     o-----o-----o  " << endl << endl;
    }




    // Lattice mesh creation
    
    latticeMesh mesh(pid);


    // Simulation handler

    timeOptions Time(pid);

    vector<uint> tlist = Time.timeList();


    // Macroscopic density

    scalarField rho( mesh, Time, "rho", IO::MUST_READ, IO::NO_WRITE );


    // Macroscopic temperature

    scalarField T( mesh, Time, "T", IO::MUST_READ, IO::NO_WRITE );


    // Macroscopic velocity

    vectorField U( mesh, Time, "U", IO::NO_READ, IO::NO_WRITE );
    

    // PDF field. Navier - Stokes equation

    pdfField f( mesh, Time, "f", IO::NO_READ, IO::NO_WRITE );


    // Macroscopic pressure

    scalarField p( mesh, Time, "p", IO::NO_READ, IO::MUST_WRITE );


    // Potential as scalar field
    
    scalarField phi( mesh, Time, "phi", IO::NO_READ, IO::MUST_WRITE );    


   
    
    
    // Navier-Stokes MRT equation

    pseudoPotEqHandler NS("Navier-Stokes", mesh, Time, f, rho, U, T);


    
    // Advance over write times
    
    for( uint i = 0 ; i < tlist.size() ; ++i ) {


	if(pid == 0) {
		
	    cout << endl <<  "Time = " << tlist[i] << endl;	       
		
	}
	

	// Update necessary fields

	rho.update(i);
	
	T.update(i);
	

	// Update potential

	NS.updatePotential(phi);


	// Compute pressure

	NS.pressure(phi, p);
	
		
    	// Write fields

	while( Time.currentTime() != tlist[i] )
	    Time.update();

	p.write();

	phi.write();
		    

    }



    






    // Print info
    if(pid == 0) {
	
    	cout << endl << "  Finished in " << Time.elapsed() << " seconds " << endl << endl;


	// Update case file

	Time.keepRegisteredFields();

	Time.updateCaseFile();
	
    }
	
    

    MPI_Finalize();

}
