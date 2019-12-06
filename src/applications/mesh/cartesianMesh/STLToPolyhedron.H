#ifndef STLTOPOLYHEDRON_H
#define STLTOPOLYHEDRON_H

#include <CGAL_datatypes.H>

Polyhedron STLToPolyhedron( const std::string& name, std::vector< std::vector<double> >& bbox );


#endif // STLTOPOLYHEDRON_H
