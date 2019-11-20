#ifndef LATTICEMESH_C_H
#define LATTICEMESH_C_H



#include "dataTypes.h"

#include "basicMesh.H"

#include "mpiInfo_C.H"
// #include <timeInfo.h>
// #include <latticeInfo.h>
// #include <EOSInfo.h>

/**
 * @file latticeMesh.h
 * @author Ezequiel O. Fogliatto
 * @date 15 Dic 2017
 * @brief Full mesh information
 */

/**
 * @brief Full mesh information
 *
 * Full mesh information
 */

typedef struct {


    // /** Time information */
    // timeInfo time;         
    

    /** Mesh information */
    basicMesh mesh;         

    /** Parallel structure */
    mpiInfo_C parallel;       

    // /** DdQq model information */
    // latticeInfo lattice;    

    // /** EOS information */
    // EOSInfo EOS;            

    uint lattice_D;
    
} latticeMesh_C;



#endif // LATTICEMESH_C_H
