#include <CGAL/IO/STL_reader.h>

#include <STLToPolyhedron.H>

#include <fstream>







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





using namespace std;

Polyhedron STLToPolyhedron( const string& name ) {


    // Load geometry

    ifstream input( name.c_str(), ios::binary);
    
    vector< CGAL::cpp11::array<double,3> > points;
    
    vector< CGAL::cpp11::array<int,3> > triangles;
    
    CGAL::read_STL(input, points, triangles);

    cout << "Reading geometry with " << points.size() << " points and "  << triangles.size() << " triangles" << endl << endl;


    
    // Construct polyhedron from triangles

    Polyhedron P;    

    for( auto tr : triangles ) {

	Build_triangle<HalfedgeDS> builder_triangle( points[tr[0]], points[tr[1]], points[tr[2]] );

	P.delegate( builder_triangle );

    }

        
    return P;

}
