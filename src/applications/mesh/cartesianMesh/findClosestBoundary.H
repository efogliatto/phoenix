#ifndef FINDCLOSESTBOUNDARY_H
#define FINDCLOSESTBOUNDARY_H

#include <vector>

#include <unordered_map>

#include <dataTypes.H>

#include <CGAL_datatypes.H>

#include <utility>

void findClosestBoundary( std::unordered_map< std::string, std::vector<uint> >& boundaries,
			  const std::vector< std::vector<uint> >& meshPoints,
			  const std::vector< std::vector<int> >& nb,
			  const std::vector< std::pair<std::string, Polyhedron> >& bdPolyMap );

#endif // FINDCLOSESTBOUNDARY_H
