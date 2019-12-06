#ifndef STLTOPOLYHEDRON_H
#define STLTOPOLYHEDRON_H

#include <Build_triangle.H>


Polyhedron STLToPolyhedron( const std::string& name, std::vector< std::vector<double> >& bbox );


#endif // STLTOPOLYHEDRON_H
