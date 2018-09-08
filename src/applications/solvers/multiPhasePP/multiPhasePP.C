#include <iostream>

#include <latticeMesh.H>

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


    
    
    latticeMesh mesh(pid);



    MPI_Finalize();

}
