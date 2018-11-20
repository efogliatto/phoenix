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
    	cout << "     o<----o---->o  Temperature gradient on non periodic boundaries" << endl;
    	cout << "     |   - | -   |  " << endl;
    	cout << "     | -   |   - |  " << endl;
    	cout << "     o-----o-----o  " << endl << endl;
    }




    // Lattice mesh creation
    
    latticeMesh mesh(pid);

    const map< string, vector<uint> >& boundaries = mesh.boundaries();



    // Open heat flux file

    if( pid == 0 ) {

	ofstream whf_file;

	whf_file.open("wallHeatFlux.dat",std::ofstream::out);

	whf_file << "#Time  ";

	for( const auto &bd : boundaries  ) {

	    whf_file << bd.first << "  ";

	}

	whf_file << endl;

	whf_file.close();

    }
    

    
    

    // Simulation handler

    timeOptions Time(pid);

    vector<uint> tlist = Time.timeList();


    // Macroscopic temperature

    scalarField T( mesh, Time, "T", IO::MUST_READ, IO::NO_WRITE );    



    // Move over time list 
    
    for( uint i = 0 ; i < tlist.size() ; ++i ) {



	// Open heat flux file

	if( pid == 0 ) {

	    ofstream whf_file;

	    whf_file.open("wallHeatFlux.dat",std::ofstream::app);

	    whf_file << tlist[i] << "  ";

	    whf_file.close();

	}
	

	
	// Read field for specific time
	
	T.update(i);


	for( const auto &bd : boundaries  ) {


	    // Compute gradient over local nodes

	    scalar grad[3] = {0,0,0};	    
	    
	    scalar g_avg[3] = {0,0,0};

	    const uint np = bd.second.size();
	    
	    
	    for( const auto j : bd.second  ) {

		T.grad( grad, j );

		for( uint k = 0 ; k < 3 ; k++ )
		    g_avg[k] += grad[k];

	    }

   

	    // Apply collective reduction

	    for( uint k = 0 ; k < 3 ; k++ ) {

		scalar globalSum = 0;

		int nelem = 0;
		

		// int nlocal = mesh.local();
    
		if( mesh.wsize() > 1 ) {

                    #ifdef DP
	    
		    MPI_Allreduce(&g_avg[k], &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                    #elif SP

		    MPI_Allreduce(&g_avg[k], &globalSum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

                    #endif
		    

		    MPI_Allreduce(&np, &nelem, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);	
	    
		}

		else {
	    
		    globalSum = g_avg[k];

		    nelem = np;

		}
		

		g_avg[k] = globalSum / nelem;
		

	    }
	    



	    


	    // Open heat flux file

	    if( pid == 0 ) {

		ofstream whf_file;

		whf_file.open("wallHeatFlux.dat",std::ofstream::app);

		whf_file << g_avg[0] << " "  << g_avg[1] << " "  << g_avg[2] << "    ";

		whf_file.close();

	    }
	    

	}




	// Open heat flux file

	if( pid == 0 ) {

	    ofstream whf_file;

	    whf_file.open("wallHeatFlux.dat",std::ofstream::app);

	    whf_file << endl;

	    whf_file.close();

	}
	
	

    }


    

    // Print info
    if(pid == 0)	
    	cout << endl << "  Finished in " << Time.elapsed() << " seconds " << endl << endl;
	
    

    MPI_Finalize();

}
