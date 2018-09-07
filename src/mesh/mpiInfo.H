#ifndef MPIINFO_H
#define MPIINFO_H

// #include <mpi.h>

#include <dataTypes.H>

#include <vector>


class mpiInfo {
    

protected:

    
    /** Process id */
    
    uint pid;   


    /** Total number of processes */
    
    uint worldSize;   


    /** Total number of ghosts */
    
    uint nghosts;   


    /** Number of local elements */
    
    uint nlocal;   


    /** Receive ghosts indices*/
    
    std::vector< std::vector<uint> > recvGhosts;   


    /** Send ghosts indices*/
    
    std::vector< std::vector<uint> > sendGhosts;           



    /** Total number of elements per patch (useful for reading/writing) */

    std::vector<uint> nodesPerPatch;

    


public:


    /** Default constructor */

    mpiInfo( const uint& id = 0 );



    /** Default destructor */

    ~mpiInfo();
    

};

#endif // MPIINFO_H
