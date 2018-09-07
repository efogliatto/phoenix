#include <iostream>

#include <latticeMesh.H>

using namespace std;


int main( int argc, char **argv ) {


    int pid, world;

    
    // Initialize mpi
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD,&pid);

    MPI_Comm_size(MPI_COMM_WORLD,&world);


    
    latticeMesh mesh(pid);



    MPI_Finalize();

}
