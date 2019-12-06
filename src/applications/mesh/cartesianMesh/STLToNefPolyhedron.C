#include <STLToNefPolyhedron.H>

#include <STLToPolyhedron.H>


using namespace std;



Nef_polyhedron STLToNefPolyhedron( const string& name, vector< vector<double> >& bbox ) {


    // Load geometry
    
    Polyhedron P = STLToPolyhedron( name, bbox );

    
    // Construct polyhedron from triangles

    Nef_polyhedron Nef(P);


        
    return Nef;

}
