#include <vectorField.H>

#include <stdio.h>


using namespace std;



/** Read field using ensight format */

const void vectorField::read() {


    
    // Open file

    if( mesh.pid() == 0 ) {

	cout << "Reading vector field " << name << endl << endl;

    }
    

    // MPI file pointer

    MPI_File file;

    MPI_File_open( MPI_COMM_WORLD,
		   ("lattice." + name + "_"+ to_string( Time.timeToIndex(Time.startTime()) )).c_str(),
		   MPI_MODE_RDONLY,
		   MPI_INFO_NULL,
		   &file );


    


    // Allocate space

    field.resize( mesh.npoints() );

    for( uint i = 0 ; i < field.size() ; i++ )
	field[i].resize( 3 );    


    float *auxField = (float*)malloc( mesh.npoints() * sizeof(float) );


    

    // Set Offset

    MPI_Offset offset = 240*sizeof(char) + sizeof(int);

    for( int i = 0 ; i < mesh.pid() ; i++ ) {

	offset += 3*mesh.npp(i) * sizeof(float);

	offset += 160*sizeof(char) + sizeof(int);

    }

    MPI_File_seek(file, offset, MPI_SEEK_SET);



    
    // Read Array

    MPI_Status st;

    for( uint j = 0 ; j < 3 ; j++ ) {

	MPI_File_read(file, auxField, mesh.npoints(), MPI_FLOAT, &st);

	for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

	    field[i][j] = (scalar)auxField[i];

	}
    
	MPI_Barrier(MPI_COMM_WORLD);

    }
    
    MPI_File_close(&file);

    free(auxField);


}


/** Default constructor */

vectorField::vectorField( const latticeMesh& m, timeOptions& t, const std::string& nm ) : latticeField(m,t,nm) {

    
    // Read values from file

    vectorField::read();
    


    // Add to time list

    Time.addVectorField(name);
    

}


/** Default destructor */

vectorField::~vectorField() {}



/** Synchronization across procceses */

const void vectorField::sync() {}




/** Write field using ensight format */

const void vectorField::write() const {


    // Open file
    
    MPI_Barrier(MPI_COMM_WORLD);
	
    string fname = "lattice." + name + "_" + to_string( Time.timeToIndex( Time.currentTime() )  );

    MPI_File file;

    MPI_File_open( MPI_COMM_WORLD,
    		   fname.c_str(),
    		   MPI_MODE_CREATE|MPI_MODE_WRONLY,
    		   MPI_INFO_NULL,
    		   &file );	

	

    // Header
	
    if( mesh.pid() == 0 ) {

    	char* msg = (char*)malloc( 80*sizeof(char) );

    	memset(msg,'\0', 80);

    	sprintf(msg, "%s", name.c_str());

    	MPI_File_write(file, msg, 80, MPI_CHAR, MPI_STATUS_IGNORE);

    	free(msg);

    }


    MPI_Barrier(MPI_COMM_WORLD);



       



    // Set Offset

    uint offset = 80*sizeof(char);

    for( int i = 0 ; i < mesh.pid() ; i++ ) {

    	offset += 3 * mesh.npp(i) * sizeof(float);

    	offset += 160*sizeof(char) + sizeof(int);

    }

    MPI_File_seek(file, offset, MPI_SEEK_SET);



    // Write "part" description

    char* msg = (char*)malloc( 80*sizeof(char) );

    memset(msg,'\0', 80);

    sprintf(msg, "part");

    MPI_File_write(file, msg, 80, MPI_CHAR, MPI_STATUS_IGNORE);
	
	
    uint pid = mesh.pid()+1;

    MPI_File_write(file, &pid, 1, MPI_INT, MPI_STATUS_IGNORE);


	
    memset(msg,'\0', 80);

    sprintf(msg, "coordinates");

    MPI_File_write(file, msg, 80, MPI_CHAR, MPI_STATUS_IGNORE);

    free(msg);
	

	
    // Write Array

    float *auxField = (float*)malloc( mesh.npoints() * sizeof(float) );

    for( uint j = 0 ; j < 3 ; j++) {

    	for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

    	    auxField[i] = (float)field[i][j];

    	}

    	MPI_File_write(file, auxField, mesh.npoints(), MPI_FLOAT, MPI_STATUS_IGNORE);

    }

    free(auxField);



    // Close file

    MPI_File_close(&file);

}
