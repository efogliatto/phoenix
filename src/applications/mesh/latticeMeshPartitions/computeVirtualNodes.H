#ifndef COMPUTEVIRTUALNODES_H
#define COMPUTEVIRTUALNODES_H

#include "../meshInclude/basicMesh.H"

#include <vector>

#include <latticeModelCreator.H>


void computeVirtualNodes( const basicMesh& mesh, std::vector< std::vector<int> >& virtualNodes, latticeModel* lbmodel );

#endif // COMPUTEVIRTUALNODES_H
