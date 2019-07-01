#include <iostream>

#include <iomanip>

#include <pseudoPotEqHandler.H>

#include <energyEqHandler.H>

#include <TEquation.H>

#include "neq.H"

#include "gammaHat.H"

#include <algebra.H>


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


    // Macroscopic velocity

    vectorField U( mesh, Time, "U", IO::MUST_READ, IO::MUST_WRITE );
    

    // PDF field. Navier - Stokes equation

    pdfField f( mesh, Time, "f", IO::MUST_READ, IO::MUST_WRITE ); 
    
    


    // Macroscopic temperature

    scalarField T( mesh, Time, "T", IO::MUST_READ, IO::MUST_WRITE );

    scalarField Tstar( mesh, Time, "T", IO::MUST_READ, IO::NO_WRITE );

    for( uint i = 0 ; i < mesh.local() ; i++ ) {

	if( mesh.latticePoint(i)[1] == 0 ) {

	    T[i] = 0.0296296;

	    Tstar[i] = 0.0296296;	    

	}

    }

    T.sync();

    Tstar.sync();




    /** Relaxation coefficients */
    
    vector<scalar> Tau;

    vector<scalar> Tau2;

    {

	dictionary dict("properties/macroProperties");

	Tau = dict.lookUp< vector<scalar> >( "Energy/LBModel/Tau" );

	Tau2 = dict.lookUp< vector<scalar> >( "Energy/LBModel/Tau" );

	for(uint k = 0 ; k < 9 ; k++) {
	    
	    Tau2[k] = 0.5 - (1.0/Tau[k]);

	}

    }



    // Non diagonal Q (0.5 - inv(Q))

    sparseScalarMatrix invQ( Tau2 );

    invQ.addElement( 0.5*(Tau[3]-1.0)/Tau[3], 3, 4);

    invQ.addElement( 0.5*(Tau[5]-1.0)/Tau[5], 5, 6);  
    
	
    


   
    // Navier-Stokes MRT equation

    pseudoPotEqHandler NS("Navier-Stokes", mesh, Time, f, rho, U, T);



    // Create eos

    EOSCreator creator;

    EOS* eos = creator.create("properties/macroProperties", "Navier-Stokes");
    

    
    
    // Advance in time. Collide, stream, update and write
    
    while( Time.update() ) {
              


	// Energy equation

	{

	    // Lattice constants

	    const vector< vector<int> >& nb = mesh.nbArray();
	    
	    const scalarMatrix& invM = mesh.lmodel()->MRTInvMatrix();

	    vector<uint> reverse = mesh.lmodel()->reverse();	    

	    const uint q = mesh.lmodel()->q();

	    vector<scalar> n_eq(q);

	    vector<scalar> n(q);
	    


	    // Predictor step

	    for( uint id = 0 ; id < mesh.local() ; id++ ) {

		if(  ( mesh.latticePoint(id)[1] > 0 )  &&  ( mesh.latticePoint(id)[1] < 300 )  ) {

		    for( uint k = 0 ; k < q ; k++ ) {

			int nbid = nb[id][k];

			if(nbid != -1)
			    n_eq[k] = neq(mesh, T, U, 1, 1, nbid, k);

			if( k == 0 )
			    n_eq[k] += 0.5*gammaHat(mesh, rho, T, U, 1, 1, Tau[3], eos, 1, id);

		    }


		    invM.matDotVec(n_eq, n);



		    Tstar[id] = 0;

		    for( uint k = 0 ; k < q ; k++ )
			Tstar[id] += n[k];

		}
		
	    }


	    Tstar.sync();



	    
	    // Corrector step

	    for( uint id = 0 ; id < mesh.local() ; id++ ) {

	    	if(  ( mesh.latticePoint(id)[1] > 0 )  &&  ( mesh.latticePoint(id)[1] < 300 )  ) {

	    	    for( uint k = 0 ; k < q ; k++ ) {

	    		int nbid = nb[id][k];

	    		int nbplus = nb[id][reverse[k]];

	    		if(  (nbid != -1)  &&  (nbplus != -1)  ) {
			
	    		    n_eq[k] = neq(mesh, Tstar, U, 1, 1, nbplus, k)
	    		        - neq(mesh, Tstar, U, 1, 1, id, k)
	    		        + neq(mesh, T, U, 1, 1, nbid, k)
	    		        - neq(mesh, T, U, 1, 1, id, k);

	    		}


	    	    }


	    	    invQ.matDotVec(n_eq, n);

	    	    invM.matDotVec(n, n_eq);

		   
		    


	    	    T[id] = 0;

	    	    for( uint k = 0 ; k < q ; k++ )
	    	    	T[id] += n_eq[k];

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
