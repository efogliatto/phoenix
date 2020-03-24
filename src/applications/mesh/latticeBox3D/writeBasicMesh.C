#include "writeBasicMesh.H"

#include <iostream>

using namespace std;

void writeBasicMesh( basicMesh& mesh ) {


    // Write files
    
    int status = system( "rm -rf lattice" );
    
    status = system( "mkdir -p lattice" );
    
    if(status) {}


    
    FILE *outFile;

    
    // Points

    outFile = fopen("lattice/points","w");   

    fprintf(outFile,"%d\n",mesh.nPoints);

    for( uint i = 0 ; i < mesh.nPoints ; i++ )
    	fprintf(outFile,"%d %d %d\n",mesh.points[i][0],mesh.points[i][1],mesh.points[i][2]);

    
    fclose(outFile);

    


    
    // Neighbours
    
    outFile = fopen("lattice/neighbours","w");

    for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

	for( uint j = 0 ; j < mesh.Q ; j++ ) {

	    fprintf(outFile,"%d ",mesh.nb[i][j]);

    	}

	fprintf(outFile,"\n");

    }

    fclose(outFile);



    


    // Boundary
    
    outFile = fopen("lattice/boundary","w");
	
    fprintf(outFile,"%d\n\n", mesh.bd.nbd);

    for( uint i = 0 ; i < mesh.bd.nbd ; i++ ) {

    	fprintf(outFile, "%s\n", mesh.bd.bdNames[i].c_str() );;
	
    	fprintf(outFile,"%d\n",mesh.bd.nbdelem[i]);

    	for( uint j = 0 ; j < mesh.bd.nbdelem[i] ; j++ ) {

    	    fprintf(outFile,"%d\n",mesh.bd.bdPoints[i][j]);

    	}

    	fprintf(outFile,"\n");
	

    }
	
	
    fclose(outFile);





    // Write VTK cells
    
    outFile = fopen("lattice/vtkCells","w");
    
    fprintf(outFile,"%d %d\n", mesh.ncells, mesh.cellType);
    
    for( uint i = 0 ; i < mesh.ncells ; i++ ) {
	
    	for( uint j = 0 ; j < mesh.cellType ; j++ ) {
	    
    	    fprintf(outFile,"%d ",mesh.vtkCells[i][j]);
	    
    	}
	
    	fprintf(outFile,"\n");
	
    }
    
    fclose(outFile);
    
}
