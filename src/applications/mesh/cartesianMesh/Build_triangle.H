#ifndef BUILD_TRIANGLE_H
#define BUILD_TRIANGLE_H


#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <CGAL/Polyhedron_3.h>


typedef CGAL::Simple_cartesian<double>     Kernel;

typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

typedef Polyhedron::HalfedgeDS             HalfedgeDS;


// A modifier creating a triangle with the incremental builder.
template <class HDS>
class Build_triangle : public CGAL::Modifier_base<HDS> {

protected:

    HalfedgeDS::Vertex::Point _p0;

    HalfedgeDS::Vertex::Point _p1;

    HalfedgeDS::Vertex::Point _p2;    
    
public:
    
    Build_triangle( CGAL::cpp11::array<double,3> point_0,
		    CGAL::cpp11::array<double,3> point_1,
		    CGAL::cpp11::array<double,3> point_2 )

	:  _p0( point_0[0], point_0[1], point_0[2] ),
	   _p1( point_1[0], point_1[1], point_1[2] ),
	   _p2( point_2[0], point_2[1], point_2[2] ) {}	


    
    void operator()( HDS& hds) {
	
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( 3, 1, 6);
        B.add_vertex( _p0 );
        B.add_vertex( _p1 );
        B.add_vertex( _p2 );
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.end_surface();
    }
};



#endif // BUILD_TRIANGLE_H
