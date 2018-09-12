#include <pdfField.H>

using namespace std;



/** Read field using ensight format */

const void pdfField::read() {


    // Open file

    if( mesh.pid() == 0 ) {

	cout << "Reading pdf field " << name << endl << endl;

    }
    


    // Allocate space

    field.resize( mesh.npoints() );

    uint q = mesh.lmodel()->q();

    for( uint i = 0 ; i < field.size() ; i++ )
	field[i].resize(q);	


    float *auxField = (float*)malloc( mesh.npoints() * sizeof(float) );


    

  
    
    for( uint k = 0 ; k < q ; k++ ) {
    

	// MPI file pointer

	MPI_File file;

	MPI_File_open( MPI_COMM_WORLD,
		       ("lattice." + name + to_string(k) + "_" + to_string( Time.timeToIndex(Time.startTime()) )).c_str(),
		       MPI_MODE_RDONLY,
		       MPI_INFO_NULL,
		       &file );



    	// Set Offset

    	MPI_Offset offset = 240*sizeof(char) + sizeof(int);

    	for( int i = 0 ; i < mesh.pid() ; i++ ) {

    	    offset += mesh.npp(i) * sizeof(float);

    	    offset += 160*sizeof(char) + sizeof(int);

    	}

    	MPI_File_seek(file, offset, MPI_SEEK_SET);



    
    	// Read Array

    	MPI_Status st;

    	MPI_File_read(file, auxField, mesh.npoints(), MPI_FLOAT, &st);

    	for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

    	    field[i][k] = (scalar)auxField[i];

    	}	
    
    	MPI_Barrier(MPI_COMM_WORLD);
    
    	MPI_File_close(&file);


    }




    free(auxField);
    

}
    

   
/** Default constructor */

pdfField::pdfField( const latticeMesh& m, timeOptions& t, const std::string& nm, const IO iopt, const IO oopt ) : latticeField(m, t, nm){


    switch(iopt) {
	

    case IO::MUST_READ:
    
	// Read values from file

	pdfField::read();

	break;


    case IO::NO_READ:

	
	// Only allocate space

	field.resize( mesh.npoints() );

	for( uint i = 0 ; i < field.size() ; i++ )
	    field[i].resize( mesh.lmodel()->q() );

	
	break;


    default:

	break;

    }
    


    

    switch (oopt) {


    case IO::MUST_WRITE:


	// Add to time list

	Time.addPdfField(name);

	break;


    default:

	break;

	
    }

}


/** Default destructor */

pdfField::~pdfField() {}



/** Synchronization across procceses */

const void pdfField::sync() {}



/** Write field using ensight format */

const void pdfField::write() const {


    float *auxField = (float*)malloc( mesh.npoints() * sizeof(float) );

    const uint q = mesh.lmodel()->q();
    

    for( uint k = 0 ; k < q ; k++ ) {


	// Open file
		
	MPI_Barrier(MPI_COMM_WORLD);

	string fname = "lattice." + name + to_string(k) + "_" + to_string( Time.timeToIndex( Time.currentTime() )  );		    

	MPI_File file;

	MPI_File_open( MPI_COMM_WORLD, fname.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file );	

	

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

	    offset += mesh.npp(i) * sizeof(float);

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

	for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

	    auxField[i] = (float)field[i][k];

	}
	    	    
	MPI_File_write(file, auxField, mesh.npoints(), MPI_FLOAT, MPI_STATUS_IGNORE);	



	// Close file

	MPI_File_close(&file);


    }


    free(auxField);
    

}
