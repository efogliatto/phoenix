#include <writeBasicMesh.H>

#include <fstream>

using namespace std;

void writeBasicMesh( const std::vector< std::vector<uint> >& points,
		     const std::vector< std::vector<int> >& nb,
		     const std::vector< std::vector<uint> >& cells,
		     const std::unordered_map< std::string, std::vector<uint> >& boundaries ) {


    // uint i,j;

    // Write files
    int status = system( "rm -rf lattice" );
    status = system( "mkdir -p lattice" );
    if(status) {}




    
    // Points

    ofstream outFile;

    outFile.open( "lattice/points" );

    outFile << points.size() << endl;

    for( uint i = 0 ; i < points.size() ; i++ )
	outFile << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;


    outFile.close();
    


    
    // Neighbours
    
    outFile.open("lattice/neighbours");

    for( uint i = 0 ; i < nb.size() ; i++ ) {

    	for( uint j = 0 ; j < nb[i].size() ; j++ ) {

    	    outFile << nb[i][j] << " ";

    	}

        outFile << endl;

    }

    outFile.close();



    


    // Boundary
    
    outFile.open("lattice/boundary");
	
    outFile << boundaries.size() << endl << endl;

    for( auto bd : boundaries ) {

        outFile << bd.first << endl;
	
    	outFile << bd.second.size() << endl;

    	for( auto pt : bd.second ) {
	    
	    outFile << pt << endl;

	}

    	outFile << endl;
	
    }
	
	
    outFile.close();





    // Write VTK cells
    
    outFile.open("lattice/vtkCells");
    
    outFile << cells.size() << " " << cells[0].size() << endl;
    						 
	
    for( auto cell : cells ) {
    
        for( auto pt : cell ) {
    
    	   outFile << pt << " ";
    
        }
    
        outFile << endl;
    
    }
    
    outFile.close();
    

    
}
