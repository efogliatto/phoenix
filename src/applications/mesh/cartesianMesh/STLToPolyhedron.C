#include <CGAL/IO/STL_reader.h>

#include <STLToPolyhedron.H>

#include <fstream>


using namespace std;



Polyhedron STLToPolyhedron( const string& name, vector< vector<double> >& bbox ) {


    // Load geometry

    ifstream input( name.c_str(), ios::binary);
    
    vector< CGAL::cpp11::array<double,3> > points;
    
    vector< CGAL::cpp11::array<int,3> > triangles;
    
    CGAL::read_STL(input, points, triangles);

    cout << "Reading geometry with " << points.size() << " points and "  << triangles.size() << " triangles" << endl << endl;



    // Out to OFF format

    ofstream output( (name + ".off").c_str() );

    output << points.size() << " " << triangles.size() << endl;

    for(auto p : points) {
	
	output << p[0] << " " << p[1] << " " << p[2] << endl;


	// Check bounding box
	
	if(p[0] < bbox[0][0])
	    bbox[0][0] = p[0];

	if(p[0] > bbox[1][0])
	    bbox[1][0] = p[0];

	if(p[1] < bbox[0][1])
	    bbox[0][1] = p[1];

	if(p[1] > bbox[1][1])
	    bbox[1][1] = p[1];

	if(p[2] < bbox[0][2])
	    bbox[0][2] = p[2];

	if(p[1] > bbox[1][2])
	    bbox[1][2] = p[2];		

    }
	

    for(auto t : triangles)
	output << 3 << " " <<  t[0] << " " << t[1] << " " << t[2] << endl;    

    
    
    // Construct polyhedron from triangles

    Polyhedron P;

    ifstream P_off( (name + ".off").c_str() );

    P_off >> P;

    // for( auto tr : triangles ) {

    // 	Build_triangle<HalfedgeDS> builder_triangle( points[tr[0]], points[tr[1]], points[tr[2]] );

    // 	P.delegate( builder_triangle );

    // }

        
    return P;

}
