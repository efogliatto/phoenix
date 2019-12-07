#include <CGAL/IO/STL_reader.h>

#include <STLToPolyhedron.H>

#include <fstream>


using namespace std;



Polyhedron STLToPolyhedron( const string& name, vector< vector<double> >& bbox ) {


    // Load geometry

    ifstream input( name.c_str() );
    
    vector< CGAL::cpp11::array<double,3> > points;
    
    vector< CGAL::cpp11::array<int,3> > triangles;
    
    CGAL::read_STL(input, points, triangles);

    cout << "Reading geometry with " << points.size() << " points and "  << triangles.size() << " triangles" << endl << endl;


    double minX(100000),
	maxX(-100000),
	minY(100000),
	maxY(-100000),
	minZ(100000),
	maxZ(-100000);
	
    

    // Out to OFF format

    ofstream output( (name + ".off").c_str() );

    output << "OFF" << endl;

    output << points.size() << " " << triangles.size() << " 0" << endl << endl;

    for(auto p : points) {
	
	output << p[0] << " " << p[1] << " " << p[2] << endl;


	// Check bounding box
	
	if(p[0] < minX)
	    minX = p[0];

	if(p[0] > maxX)
	    maxX = p[0];

	if(p[1] < minY)
	    minY = p[1];

	if(p[1] > maxY)
	    maxY = p[1];

	if(p[2] < minZ)
	    minZ = p[2];

	if(p[2] > maxZ)
	    maxZ = p[2];		

    }


    bbox = { {minX, minY, minZ}, {maxX, maxY, maxZ} };
	

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
