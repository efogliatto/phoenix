#include <cellsInsidePolyhedron.H>

#include <CGAL/AABB_tree.h>

#include <CGAL/AABB_traits.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/algorithm.h>

#include <CGAL/Side_of_triangle_mesh.h>


typedef Kernel::Point_3 Point;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;

typedef CGAL::AABB_traits<Kernel, Primitive> Traits;

typedef CGAL::AABB_tree<Traits> Tree;

typedef CGAL::Side_of_triangle_mesh<Polyhedron, Kernel> Point_inside;



using namespace std;




void cellsInsidePolyhedron( vector<bool>& isInside, const Polyhedron& P, const vector< vector<uint> >& points, const vector< vector<uint> >& cells ) {


    // Construct AABB tree with a KdTree

    Tree tree( faces(P).first, faces(P).second, P );
    
    tree.accelerate_distance_queries();


    // Initialize the point-in-polyhedron tester
    
    Point_inside inside_tester(tree);


    // Allocate space for isInside

    isInside.resize( cells.size() );


    
    // For each cell, check if center of mass is inside Polyhedron

    for( uint i = 0 ; i < cells.size() ; i++) {


	double com[3] = {0,0,0};

	for( uint j = 0 ; j < cells[i].size() ; j++ ) {

	    com[0] += points[ cells[i][j] ][0];
	    
	    com[1] += points[ cells[i][j] ][1];

	    com[2] += points[ cells[i][j] ][2];
	    
	}


	for( uint j = 0 ; j < 3 ; j++ )
	    com[j] = com[j] / cells[i].size();




	// Create point and check if is inside poly

	Point query( com[0], com[1], com[2] );
	
	isInside[i] = ( inside_tester(query) == CGAL::ON_BOUNDED_SIDE );
	

    }

    
    

}
