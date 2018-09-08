#include <scalarField.H>

using namespace std;


/** Read field using ensight format */

const void scalarField::read() {

    if( mesh.pid() == 0 ) {

	cout << "Reading scalar field " << name << endl << endl;

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



    // Set Offset

    MPI_Offset offset = 240*sizeof(char) + sizeof(int);

    for( int i = 0 ; i < mesh.pid() ; i++ ) {

    	offset += mesh.npp(i) * sizeof(float);

    	offset += 160*sizeof(char) + sizeof(int);

    }

    MPI_File_seek(file, offset, MPI_SEEK_SET);



    
    // Read Array

    float* auxField = (float*)malloc( mesh.npoints() * sizeof(float) );
	
    MPI_Status st;

    MPI_File_read(file, auxField, mesh.npoints(), MPI_FLOAT, &st);

    for( uint i = 0 ; i < mesh.npoints() ; i++) {

    	field[i] = (scalar)auxField[i];
	    
    }	
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_File_close(&file);

    free(auxField);
    
    

}









/** Default constructor */

scalarField::scalarField( const latticeMesh& m, const timeOptions& t, const string& nm ) : latticeField(m,t,nm) {

    
    // Read values from file

    scalarField::read();
    

}


/** Default destructor */

scalarField::~scalarField() {}



// Acces members. Forced interface

/** Synchronization across procceses */

const void scalarField::sync() {}


/** Write field using ensight format */

const void scalarField::write() const {}

