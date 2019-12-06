#ifndef STLTOPOLYHEDRON_H
#define STLTOPOLYHEDRON_H

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <CGAL/Polyhedron_3.h>


typedef CGAL::Simple_cartesian<double>     Kernel;

typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

typedef Polyhedron::HalfedgeDS             HalfedgeDS;


Polyhedron STLToPolyhedron( const std::string& name );


#endif // STLTOPOLYHEDRON_H
