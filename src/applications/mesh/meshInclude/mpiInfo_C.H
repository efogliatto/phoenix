#ifndef MPIINFO_C_H
#define MPIINFO_C_H

/**< MPI buffer size */
#define MPI_BUFF_SIZE 500  

/**
 * @file mpiInfo.h
 * @author Ezequiel O. Fogliatto
 * @date 15 Dic 2017
 * @brief Struct containing mpi information
 */

/**
 * @brief Struct containing mpi information
 *
 * This struct contains information related to mpi communications
 */

#include <dataTypes.H>

#include <vector>

typedef struct {

    
    /**< Process id */
    
    uint pid;   


    /**< Total number of processes */
    
    uint worldSize;   


    /**< Total number of ghosts */
    
    uint nghosts;   


    /**< Number of local elements */
    
    uint nlocal;   


    /**< Receive ghosts indices*/
    
    std::vector< std::vector<uint> > recvGhosts;   


    /**< Send ghosts indices*/
    
    std::vector< std::vector<uint> > sendGhosts;           


    /**< Number of receive ghosts */
    
    std::vector<uint> nrg;   


    /**< Number of send ghosts */
    
    std::vector<uint> nsg;


    /** Total number of elements per patch (useful for reading/writing) */

    std::vector<uint> nodesPerPatch;
    

    
    // Buffers

    /** Scalar send buffer */
    
    std::vector< std::vector<scalar> > sbuf;   


    /** Scalar receive buffer */
    
    std::vector< std::vector<scalar> > rbuf;   

    
    /** Vector send buffer */
    
    std::vector< std::vector<scalar> > vsbuf;   

    
    /** Vector receive buffer */
    
    std::vector< std::vector<scalar> > vrbuf;



    /** Grid node level */

    std::vector<uint> level;
    
    
} mpiInfo_C;

#endif // MPIINFO_C_H
