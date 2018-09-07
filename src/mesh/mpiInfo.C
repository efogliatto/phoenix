#include <mpiInfo.H>

#include <fstream>

#include <sstream>

#include <string>


using namespace std;


/** Default constructor */

mpiInfo::mpiInfo( const uint& id ) : pid(id) {


    // // Read recv ghosts
    
    // ostringstream fname;
    // fname << "processor" << pid << "/lattice/recvGhosts";
    // ifstream inFile( fname.c_str() );
    

    // sprintf(command,"processor%d/lattice/recvGhosts", pid);

    // inFile = fopen( command, "r" );

  
    // status = fscanf(inFile, "%d\n", &mesh.parallel.worldSize);

    // mesh.parallel.nrg        = (uint*)malloc( mesh.parallel.worldSize * sizeof(uint) );

    // mesh.parallel.recvGhosts = (uint**)malloc( mesh.parallel.worldSize * sizeof(uint*) );
    
   
    
    // for( i = 0 ; i < mesh.parallel.worldSize ; i++ ) {

    // 	uint auxPid;
	
    // 	status = fscanf(inFile, "%d", &auxPid );


    // 	if( auxPid == i ) {
	

    // 	    status = fscanf(inFile, "%d", &mesh.parallel.nrg[i] );
	
    // 	    mesh.parallel.recvGhosts[i] = (uint*)malloc( mesh.parallel.nrg[i] * sizeof(uint) );

	
    // 	    for( j = 0 ; j < mesh.parallel.nrg[i] ; j++ ) {

    // 		status = fscanf( inFile, "%d", &mesh.parallel.recvGhosts[i][j] );

    // 	    }

	    
    // 	}

    // }

    
    // fclose(inFile);

    

}



/** Default destructor */

mpiInfo::~mpiInfo() {}
