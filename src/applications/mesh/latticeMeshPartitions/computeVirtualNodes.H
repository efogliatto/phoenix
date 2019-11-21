#ifndef COMPUTEVIRTUALNODES_H
#define COMPUTEVIRTUALNODES_H


#include "../meshInclude/basicMesh.H"

#include <vector>

#include <latticeModelCreator.H>


/**
 * @file computeVirtualNodes.H
 * @author Ezequiel O. Fogliatto
 * @date 21 Nov 2019
 * @brief Compute virtual nodes for boundary nodes
 */

/**
 * Computes virtual nodes for boundary nodes, and finds connection for missing neighbours.
 * For every missing velocity, find closest node on boundary and it neighbour along the same direction
 * @p virtualNodes Node array (not allocated). Each element contains 4 indices: base node, velocity, closest node on boundary (closest to virtual), neighbour to closest node
 * @return Mesh structure
 */

void computeVirtualNodes( const basicMesh& mesh, std::vector< std::vector<int> >& virtualNodes, latticeModel* lbmodel );

#endif // COMPUTEVIRTUALNODES_H
