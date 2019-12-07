#ifndef CELLSINSIDEPOLYHEDRON_H
#define CELLSINSIDEPOLYHEDRON_H

#include <CGAL_datatypes.H>

void cellsInsidePolyhedron( std::vector<bool>& isInside, const Polyhedron& P, const std::vector< std::vector<uint> >& points, const std::vector< std::vector<uint> >& cells );

#endif // CELLSINSIDEPOLYHEDRON_H
