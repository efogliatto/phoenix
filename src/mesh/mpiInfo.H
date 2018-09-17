#ifndef MPIINFO_H
#define MPIINFO_H

#include <mpi.h>

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


    /** Processor id */

    const uint id() const {  return pid; }

    
    /** Total number of processors */

    const uint wsize() const {  return worldSize; }
   

    /** Number of local nodes */

    const uint local() const {  return nlocal; }


    /** Number of ghost nodes */

    const uint ghosts() const {  return nghosts; }


    /** Nodes per patch */

    const uint npp( const uint& id ) const { return nodesPerPatch[id]; }


    /** Reference to recvGhosts */

    const std::vector< std::vector<uint> >& recvg() const { return recvGhosts; }


    /** Reference to sendGhosts */

    const std::vector< std::vector<uint> >& sendg() const { return sendGhosts; }    

};

#endif // MPIINFO_H
