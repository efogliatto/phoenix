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

    // scalarField rho( mesh, Time, "rho", IO::MUST_READ, IO::NO_WRITE );


    // Macroscopic temperature

    scalarField T( mesh, Time, "T", IO::MUST_READ, IO::NO_WRITE );

    scalarField oldT( mesh, Time, "T", IO::NO_READ, IO::NO_WRITE );

    scalarField residual( mesh, Time, "residual", IO::NO_READ, IO::MUST_WRITE );
    


    // Macroscopic velocity

    vectorField U( mesh, Time, "U", IO::NO_READ, IO::NO_WRITE );


    // // Temperature laplacian

    // scalarField lapT( mesh, Time, "laplacian_T", IO::NO_READ, IO::MUST_WRITE );


    // Temperature gradient

    vectorField gradT( mesh, Time, "gradient_T", IO::NO_READ, IO::MUST_WRITE );


   
     

    
    // Advance over write times
    
    for( uint i = 0 ; i < tlist.size() ; ++i ) {


	if(pid == 0) {
		
	    cout << endl <<  "Time = " << tlist[i] << endl;	       
		
	}
	

	// Update necessary fields
	
	T.update(i);

	U.update(i);


	
	// Compute v * T

	for( uint id = 0 ; id < mesh.local() ; id++ ) {

	    for( uint j = 0 ; j < 3 ; j++ )
		gradT[id][j] = U.at(id)[j] * T.at(id);

	}

	gradT.sync();
	

	
	// Compute residuals

	for( uint id = 0 ; id < mesh.local() ; id++ ) {

	    residual[id] = T.at(id) - oldT.at(id)
		         + gradT.div(id)
		         - 0.5 * T.laplacian(id) / 3;

	}

	residual.sync();


	for( uint id = 0 ; id < mesh.local() ; id++ )
	    oldT[id] = T.at(id);

	    

	
		
    	// Write fields

	while( Time.currentTime() != tlist[i] )
	    Time.update();

	residual.write();
	


		    

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
