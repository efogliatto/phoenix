#ifndef WRITEBASICMESH_H
#define WRITEBASICMESH_H

#include <vector>

#include <map>

#include <string>

void writeBasicMesh( const std::vector< std::vector<uint> >& points,
		     const std::vector< std::vector<int> >& nb,
		     const std::vector< std::vector<uint> >& cells,
		     const std::map< std::string, std::vector<uint> >& boundaries );

#endif // WRITEBASICMESH_H
