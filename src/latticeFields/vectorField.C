#include <vectorField.H>

#include <stdio.h>

#include <cmath>


using namespace std;




/** Default constructor */

vectorField::vectorField( const latticeMesh& m, timeOptions& t, const std::string& nm, const IO iopt, const IO oopt ) : latticeField(m,t,nm) {


    switch(iopt) {
	

    case IO::MUST_READ:    

    
	// Read values from file

	vectorField::read();

	break;



    case IO::NO_READ:

	
	// Only allocate space

	field.resize( mesh.npoints() );

	for( uint i = 0 ; i < field.size() ; i++ )
	    field[i].resize( 3 );   	
	
	break;


    default:

	break;



    }






    switch (oopt) {


    case IO::MUST_WRITE:	
    

	// Add to time list

	Time.addVectorField(name);

	break;


    default:

	break;


    }
    



    // Resize mpi send buffers

    const vector<vector<uint>>& sendGhosts = mesh.mpi().sendg();
    
    sbuf = (scalar**)malloc( mesh.wsize() * sizeof(scalar*) );
    
    for( int i = 0 ; i < mesh.wsize() ; i++ )
    	sbuf[i] = (scalar*)malloc( sendGhosts[i].size() * 3 * sizeof(scalar) );


    

    // Resize mpi recv buffers
    
    const vector<vector<uint>>& recvGhosts = mesh.mpi().recvg();

    rbuf = (scalar**)malloc( mesh.wsize() * sizeof(scalar*) );
    
    for( int i = 0 ; i < mesh.wsize() ; i++ ) 
    	rbuf[i] = (scalar*)malloc( recvGhosts[i].size() * 3 * sizeof(scalar) );

    
}



/** Default destructor */

vectorField::~vectorField() {


    if( mesh.wsize() > 1 ) {


	for( int i = 0 ; i < mesh.wsize() ; i++ ) {
	    
	    free(sbuf[i]);

	    free(rbuf[i]);

	}

	free(sbuf);

	free(rbuf);

    }
    

}





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







/** Write field using ensight format */

const void vectorField::write() const {


    // First check for array sanity   

    for( uint i = 0 ; i < mesh.npoints() ; i++) {

	for( uint j = 0 ; j < 3 ; j++) {

	    if( isnan(field[i][j]) ) {

		cout << " [ERROR] Floating point exception. NaN solution" << endl;

		exit(1);

	    }

	}

    }


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






/** Synchronization across procceses */

const void vectorField::sync() {


    if( mesh.wsize() > 1 ) {


	// Reference to mpi info

	const vector<vector<uint>>& sendGhosts = mesh.mpi().sendg();

	const vector<vector<uint>>& recvGhosts = mesh.mpi().recvg();	

	
	
	// Move over send ghosts. Copy data to buffer

	for( int pid = 0 ; pid < mesh.wsize() ; pid++ ) {


	    for( uint i = 0 ; i < sendGhosts[pid].size() ; i++) {


	    	for( uint k = 0 ; k < 3 ; k++ ) {


	    	    // Position in sbuf

	    	    uint id = i*3 + k;

		    
	    	    // Copy field info
		    
	    	    sbuf[pid][id] = field[ sendGhosts[pid][i] ][ k ];
		    

	    	}

	    }
	  

	}




    	// Send information
    	
    	MPI_Status status[100];

	MPI_Request request[100];

	int nreq = 0;

	
	// Move over send ghosts. Send data

	for( int pid = 0 ; pid < mesh.wsize() ; pid++ ) {

	    if( sendGhosts[pid].size() > 0 ) {

		#ifdef DP
		
		MPI_Isend (&sbuf[pid][0], sendGhosts[pid].size() * 3, MPI_DOUBLE, pid, mesh.pid(), MPI_COMM_WORLD, &request[nreq]);

		#elif SP

		MPI_Isend (&sbuf[pid][0], sendGhosts[pid].size() * 3, MPI_FLOAT, pid, mesh.pid(), MPI_COMM_WORLD, &request[nreq]);

		#endif

		nreq++;

	    }
	}



	// Move over recv ghosts. Receive data

	for( int pid = 0 ; pid < mesh.wsize() ; pid++ ) {

	    if( recvGhosts[pid].size() > 0 ) {

		#ifdef DP
		
		MPI_Irecv (rbuf[pid], recvGhosts[pid].size() * 3, MPI_DOUBLE, pid, pid, MPI_COMM_WORLD, &request[nreq]);

		#elif SP

		MPI_Irecv (rbuf[pid], recvGhosts[pid].size() * 3, MPI_DOUBLE, pid, pid, MPI_COMM_WORLD, &request[nreq]);

		#endif

		nreq++;

	    }
	}


	
	// Wait for everyone

	MPI_Waitall(nreq, request, status);




	

	// Copy new data back to ghost nodes

	for( int pid = 0 ; pid < mesh.wsize() ; pid++ ) {
	   

	    for( uint i = 0 ; i < recvGhosts[pid].size() ; i++) {

		for( uint k = 0 ; k < 3 ; k++ ) {


		    // Position in sbuf

		    uint id = i*3 + k;

		    
		    // Copy field info
		    
		    field[ recvGhosts[pid][i] ][ k ] = rbuf[pid][id];
		    

		}

	    }		

	    
	}	
	


    }
    

}




/** Start sync */

const void vectorField::startSync() {}


/** End sync */

const void vectorField::endSync() {}



/** Vector divergence at node */

const scalar vectorField::div( const uint& id ) const {


    // Constants

    int a, b;

    scalar d(0);

    const vector< vector<int> >& nb = mesh.nbArray();
        

    // D2Q9 model

    if( mesh.lmodel()->name() == "D2Q9" ) {

    
    	// X - derivative

    	a = nb[id][3];
	
    	b = nb[id][1];

	
    	if(  (a != -1)  &&  (b != -1)  ) {
    
    	    d += 0.5 * (field[a][0] - field[b][0]);

    	}

    	else {

    	    if(  (a == -1)  &&  (b != -1)  ) {
    
    		d += (field[id][0] - field[b][0]);

    	    }

    	    else {

    		if(  (a != -1)  &&  (b == -1)  ) {
    
    		    d += (field[a][0] - field[id][0]);

    		}

    	    }

    	}




    	// Y - derivative

    	a = nb[id][4];
	
    	b = nb[id][2];

	
    	if(  (a != -1)  &&  (b != -1)  ) {
    
    	    d += 0.5 * (field[a][1] - field[b][1]);

    	}

    	else {

    	    if(  (a == -1)  &&  (b != -1)  ) {
    
    		d += (field[id][1] - field[b][1]);

    	    }

    	    else {

    		if(  (a != -1)  &&  (b == -1)  ) {
    
    		    d += (field[a][1] - field[id][1]);

    		}

    	    }

    	} 
	

	

    }

    

    else {

	if( mesh.lmodel()->name() == "D3Q15") {


	    // X - derivative

	    a = nb[id][2];
	
	    b = nb[id][1];

	
	    if(  (a != -1)  &&  (b != -1)  ) {
    
		d += 0.5 * (field[a][0] - field[b][0]);

	    }

	    else {

		if(  (a == -1)  &&  (b != -1)  ) {
    
		    d += (field[id][0] - field[b][0]);

		}

		else {

		    if(  (a != -1)  &&  (b == -1)  ) {
    
			d += (field[a][0] - field[id][0]);

		    }

		}

	    }




	    // Y - derivative

	    a = nb[id][4];
	
	    b = nb[id][3];

	
	    if(  (a != -1)  &&  (b != -1)  ) {
    
		d += 0.5 * (field[a][1] - field[b][1]);

	    }

	    else {

		if(  (a == -1)  &&  (b != -1)  ) {
    
		    d += (field[id][1] - field[b][1]);

		}

		else {

		    if(  (a != -1)  &&  (b == -1)  ) {
    
			d += (field[a][1] - field[id][1]);

		    }

		}

	    }   




	    // Z - derivative

	    a = nb[id][6];
	
	    b = nb[id][5];

	
	    if(  (a != -1)  &&  (b != -1)  ) {
    
		d += 0.5 * (field[a][2] - field[b][2]);

	    }

	    else {

		if(  (a == -1)  &&  (b != -1)  ) {
    
		    d += (field[id][2] - field[b][2]);

		}

		else {

		    if(  (a != -1)  &&  (b == -1)  ) {
    
			d += (field[a][2] - field[id][2]);

		    }

		}

	    }   


	}

    }


	
       

    return d;

}
