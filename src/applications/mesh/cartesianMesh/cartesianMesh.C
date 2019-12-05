/*

  latticeMeshPartition

  Mesh subdivision for parallel processing

 */


#include <iostream>

#include <latticeModelCreator.H>

#include <dictionary.H>


#include <CGAL/IO/STL_reader.h>

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef Kernel::Point_2 Point_2;
// typedef Kernel::Segment_2 Segment_2;



using namespace std;


int main(int argc, char** argv) {



    cout << "                    " << endl;
    cout << "     o-----o-----o  " << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     |   - | -   |        cartesianMesh" << endl;
    cout << "     o<----o---->o  " << endl;
    cout << "     |   - | -   |  Cartesian Mesh Generator" << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     o-----o-----o  " << endl << endl;



    
    // Read model name. Use D2Q9 as default

    dictionary ldict("properties/latticeProperties");
    
    latticeModelCreator lbm;
	
    latticeModel* lbmodel = lbm.create(ldict.lookUpOrDefault<string>("LBModel","D2Q9"));



    // Load main geometry

    std::ifstream input(argv[1], std::ios::binary);
    std::vector< CGAL::cpp11::array<double,3> > points;
    std::vector< CGAL::cpp11::array<int,3> > triangles;

   CGAL::read_STL(input, points, triangles);

   std::cout << points.size() << "  "  << triangles.size() << std::endl;

    


    
    



    cout << "Finished meshing" << endl << endl;

   
    
    return 0;

}
