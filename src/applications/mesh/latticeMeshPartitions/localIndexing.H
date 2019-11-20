#ifndef LOCALINDEXING_H
#define LOCALINDEXING_H

#include <vector>

#include "../meshInclude/basicMesh.H"


void localIndexing( basicMesh& mesh, std::vector< std::vector<int> >& local, std::vector< std::vector<int> >& nGhosts, std::vector<uint>& owner, uint np );


#endif // LOCALINDEXING_H
