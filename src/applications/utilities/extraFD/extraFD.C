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
    	cout << "     o<----o---->o  Additional calculations" << endl;
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

    scalarField lapRho( mesh, Time, "lapRho", IO::NO_READ, IO::MUST_WRITE );


    // Apparent contact angle

    scalarField contactAngle( mesh, Time, "contactAngle", IO::NO_READ, IO::MUST_WRITE );

    for(uint i = 0 ; i < mesh.local() ; i++)
	contactAngle[i] = 180;

    contactAngle.sync();
    


    // // Macroscopic temperature

    // scalarField T( mesh, Time, "T", IO::MUST_READ, IO::NO_WRITE );


    // // Macroscopic velocity

    // vectorField U( mesh, Time, "U", IO::NO_READ, IO::NO_WRITE );
    

    // // PDF field. Navier - Stokes equation

    // pdfField f( mesh, Time, "f", IO::NO_READ, IO::NO_WRITE );


    // // Macroscopic pressure

    // scalarField p( mesh, Time, "p", IO::NO_READ, IO::MUST_WRITE );


    // // Potential as scalar field
    
    // scalarField phi( mesh, Time, "phi", IO::NO_READ, IO::NO_WRITE );    


   
    
    
    // // Navier-Stokes MRT equation

    // pseudoPotEqHandler NS("Navier-Stokes", mesh, Time, f, rho, U, T);


    if(pid == 0)
	cout << "Computing density laplacian and apparent contact angle" << endl;
    
    
    // Advance over write times
    
    for( uint i = 0 ; i < tlist.size() ; ++i ) {


	if(pid == 0) {
		
	    cout << endl <<  "Time = " << tlist[i] << endl;	       
		
	}
	

	// Update necessary fields

	rho.update(i);


	for(uint id = 0; id < mesh.local() ; id++)
	    lapRho[id] = rho.laplacian(id);



	// Density laplacian
	
	lapRho.sync();



	// Apparent contact angle

	for( uint i = 0 ; i < mesh.local() ; i++ ) {

	    if( mesh.latticePoint(i)[2] == 0 ) {

		scalar apangle(180);
		
		scalar gradRho[3] = {0,0,0};

		rho.cartesianGradient(gradRho, i);
				
		scalar gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] );


		// Compute aparent angle

		if(gmag != 0) {
				
		    apangle = -gradRho[2]  /  gmag;

		    apangle = M_PI/2 - atan(apangle);

		    apangle = apangle * 180 / M_PI;
			
		}

		contactAngle[i] = apangle;

	    }

	}


	
	
	
	// T.update(i);
	

	// // Update potential

	// NS.updatePotential(phi);


	// // Compute pressure

	// NS.pressure(phi, p);
	
		
    	// Write fields

	while( Time.currentTime() != tlist[i] )
	    Time.update();

	lapRho.write();

	contactAngle.write();
		    

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
