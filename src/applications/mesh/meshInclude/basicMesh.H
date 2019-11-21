#ifndef BASICMESH_H
#define BASICMESH_H

/**
 * @file basicMesh.H
 * @author Ezequiel O. Fogliatto
 * @date 15 Dic 2017
 * @brief Mesh information
 */

#include "dataTypes.h"

#include <vector>

#include <string>

/**
 * @brief Boundary information
 *
 * Information about boundary elements
 */

typedef struct {
    
    uint nbd;         /**< Number of boundary types */
    
    std::vector<uint> nbdelem;    /**< Total number of elements per boundary type */
    
    std::vector<std::string> bdNames;   /**< Boundary names */
    
    std::vector< std::vector<uint> > bdPoints;  /**< Elements in boundary */

} basicBoundary;



/**
 * @brief Basic lattice information
 *
 * Lattice information
 */

typedef struct {

    
    uint nPoints;   /**< Total number of points per patch */
    
    std::vector< std::vector<int> > points;   /**< Points coordinates array */
    

    
    uint Q;   /**< Total number of neighbours */
    
    char lbm[100];   /**< LB lattice model (DdQq type) */
    
    std::vector< std::vector<int> > nb;   /**< Neighbour indices array */


    
    basicBoundary bd;   /**< Boundary information */

    
    uint ncells;   /**< Number of VTK cells per patch */
    
    uint cellType;   /**< VTK cell type */
    
    std::vector< std::vector<int> > vtkCells;   /**< VTK cell indices array */


    /** Array with cell conectivity */

    std::vector< std::vector<uint> > nodeToCells;


  /** Flag to check if node is on boundary */
  
    uint* isOnBnd;

} basicMesh;





#endif // BASICMESH_H
