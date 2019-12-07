#ifndef CGAL_DATATYPES_H
#define CGAL_DATATYPES_H

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/Nef_polyhedron_3.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
// typedef CGAL::Simple_cartesian<double>     Kernel;

typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

typedef Polyhedron::HalfedgeDS             HalfedgeDS;

typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron;


#endif // CGAL_DATATYPES_H
