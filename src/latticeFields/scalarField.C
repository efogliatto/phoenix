#include <scalarField.H>

#include <stdio.h>

#include <cmath>


using namespace std;



/** Default constructor */

scalarField::scalarField( const latticeMesh& m, timeOptions& t, const string& nm, const IO iopt, const IO oopt ) : latticeField(m,t,nm) {


    switch(iopt) {
	

    case IO::MUST_READ:

	
	// Read values from file

	scalarField::read();

	break;


    case IO::NO_READ:

	
	// Only allocate space
        
	field.resize( mesh.npoints() );

	break;


	
    default:

	break;

    }


    

    switch (oopt) {


    case IO::MUST_WRITE:
	
	// Add to time list

	Time.addScalarField(name);

	break;

	
    default:

	break;

    }
    


    if( mesh.wsize() > 1 ) {
	

	// Resize mpi send buffers

	const vector<vector<uint>>& sendGhosts = mesh.mpi().sendg();
    
	sbuf = new scalar* [mesh.wsize()];
    
	for( int i = 0 ; i < mesh.wsize() ; i++ )
	    sbuf[i] = new scalar [sendGhosts[i].size()];


    

	// Resize mpi recv buffers
    
	const vector<vector<uint>>& recvGhosts = mesh.mpi().recvg();

	rbuf = new scalar* [mesh.wsize()];
    
	for( int i = 0 ; i < mesh.wsize() ; i++ ) 
    	rbuf[i] = new scalar [recvGhosts[i].size()];

	
    }
    
}




/** Default destructor */

scalarField::~scalarField() {

    if( mesh.wsize() > 1 ) {

	for( int i = 0 ; i < mesh.wsize() ; i++ ) {
	    
	    delete [] sbuf[i];

	    delete [] rbuf[i];

	}

	delete [] sbuf;

	delete [] rbuf;

    }

}






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




// Acces members. Forced interface

/** Synchronization across procceses */

const void scalarField::sync() {

    startSync();

    endSync();

}




/** Start sync */

const void scalarField::startSync() {


    if( mesh.wsize() > 1 ) {


	// Reference to mpi info

	const vector<vector<uint>>& sendGhosts = mesh.mpi().sendg();

	const vector<vector<uint>>& recvGhosts = mesh.mpi().recvg();	


	// Move over send ghosts. Copy data to buffer

	for( int pid = 0 ; pid < mesh.wsize() ; pid++ ) {


	    for( uint i = 0 ; i < sendGhosts[pid].size() ; i++) {

		    
		// Copy field info
		    
		sbuf[pid][i] = field[ sendGhosts[pid][i] ];


	    }
	  

	}



    	// Send information    	

	_nreq = 0;

	
	// Move over send ghosts. Send data

	for( int pid = 0 ; pid < mesh.wsize() ; pid++ ) {

	    if( sendGhosts[pid].size() > 0 ) {

		#ifdef DP
		
		MPI_Isend (&sbuf[pid][0], sendGhosts[pid].size(), MPI_DOUBLE, pid, mesh.pid(), MPI_COMM_WORLD, &_request[_nreq]);

		#elif SP

		MPI_Isend (&sbuf[pid][0], sendGhosts[pid].size(), MPI_FLOAT, pid, mesh.pid(), MPI_COMM_WORLD, &_request[_nreq]);

		#endif

		_nreq = _nreq + 1;

	    }
	}



	// Move over recv ghosts. Receive data

	for( int pid = 0 ; pid < mesh.wsize() ; pid++ ) {

	    if( recvGhosts[pid].size() > 0 ) {

		#ifdef DP
		
		MPI_Irecv (rbuf[pid], recvGhosts[pid].size(), MPI_DOUBLE, pid, pid, MPI_COMM_WORLD, &_request[_nreq]);

		#elif SP

		MPI_Irecv (rbuf[pid], recvGhosts[pid].size(), MPI_DOUBLE, pid, pid, MPI_COMM_WORLD, &_request[_nreq]);

		#endif

		_nreq = _nreq + 1;

	    }
	}
	
	

    }
    

}






/** End sync */

const void scalarField::endSync() {

    
    if( mesh.wsize() > 1 ) {


	// Mpi references
	
	const vector<vector<uint>>& recvGhosts = mesh.mpi().recvg();

	
    	// Wait for everyone

    	MPI_Waitall(_nreq, _request, _status);

	

    	// Copy new data back to ghost nodes
	
    	for( int pid = 0 ; pid < mesh.wsize() ; pid++ ) {	   

    	    for( uint i = 0 ; i < recvGhosts[pid].size() ; i++)	   		    
		field[ recvGhosts[pid][i] ] = rbuf[pid][i];
		    
	    
    	}	


    }
    

}






/** Write field using ensight format */

const void scalarField::write() const {



    // First check for array sanity   

    for( uint i = 0 ; i < mesh.npoints() ; i++) {

	if( std::isnan(field[i]) ) {

	    cout << " [ERROR] Floating point exception. NaN solution" << endl;

	    exit(1);

	}


    }
    
	
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

    	offset += mesh.npp(i) * sizeof(float);

    	offset += 160*sizeof(char) + sizeof(int);

    }

    MPI_File_seek(file, offset, MPI_SEEK_SET);


    // Write "part" description

    char* msg = (char*)malloc( 80*sizeof(char) );

    memset(msg,'\0', 80);	

    sprintf(msg, "part");

    MPI_File_write(file, msg, 80, MPI_CHAR, MPI_STATUS_IGNORE);
	
	
    int pid = mesh.pid()+1;

    MPI_File_write(file, &pid, 1, MPI_INT, MPI_STATUS_IGNORE);


	
    memset(msg,'\0', 80);	

    sprintf(msg, "coordinates");

    MPI_File_write(file, msg, 80, MPI_CHAR, MPI_STATUS_IGNORE);

    free(msg);
	

	
    // Write Array

    float *auxField = (float*)malloc( mesh.npoints() * sizeof(float) );

    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

    	auxField[i] = (float)field[i];

    }
	
    MPI_File_write(file, auxField, mesh.npoints(), MPI_FLOAT, MPI_STATUS_IGNORE);

    free(auxField);



    // Close file

    MPI_File_close(&file);
    

}





/** Gradient at node */

const void scalarField::grad(scalar g[3], const uint& id, const bool inverse) const {


    // Lattice constants

    const vector< vector<int> >& nb = mesh.nbArray();

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();

    const vector<uint>& reverse = mesh.lmodel()->reverse();

    const vector<scalar>& omega = mesh.lmodel()->omega();

    const scalar cs2 = mesh.lmodel()->cs2();

    const scalar q = mesh.lmodel()->q();


    
    // Initialize gradient

    for( uint j = 0 ; j < 3 ; j++ )
	g[j] = 0;



    // Move over velocities
    
    for( uint k = 1 ; k < q ; k++ ) {

	int nbid = nb[id][k];

	if( nbid != -1 ) {

	    for( uint j = 0 ; j < 3 ; j++ )
		inverse  ?  g[j] += omega[k] * vel[reverse[k]][j] / (field[nbid] * cs2 )  :  g[j] += omega[k] * vel[reverse[k]][j] * field[nbid] / cs2;

	}


	// If neighbour does not exist, use virtual node

	else {

	    int first( mesh.vnode(id, reverse[k]) ),
		second( mesh.vnode(id, reverse[k], false) );

	    scalar fval(0);

	    if( second != -1 ) {

		inverse ? fval = 2/field[first] - 1/field[second] : fval = 2*field[first] - field[second];

	    }

	    else {

		inverse ? fval = 1/field[first] : fval = field[first];

	    }


	    // Add virtual node contribution

	    for( uint j = 0 ; j < 3 ; j++ )
		g[j] += omega[k] * vel[reverse[k]][j] * fval / cs2;


	    

	    // int otherNb = nb[id][reverse[k]];

	    // if( otherNb != -1 ) {

	    // 	for( uint j = 0 ; j < 3 ; j++ )
	    // 	    inverse  ?  g[j] += omega[k] * vel[reverse[k]][j] / ((2.0*field[id] - field[otherNb]) * cs2 )  :  g[j] += omega[k] * vel[reverse[k]][j] * (2.0*field[id]-field[otherNb]) / cs2;

	    // }

	    // else {

	    // 	for( uint j = 0 ; j < 3 ; j++ )
	    // 	    inverse  ?  g[j] += omega[k] * vel[reverse[k]][j] / (field[id] * cs2 )  :  g[j] += omega[k] * vel[reverse[k]][j] * field[id] / cs2;

	    // }

	    

	}

	

    }
    

}












/** Laplacian at node */

const scalar scalarField::laplacian(const uint& id, const bool inverse) const {


    // Lattice constants

    const vector< vector<int> >& nb = mesh.nbArray();

    const vector<uint>& reverse = mesh.lmodel()->reverse();

    const vector<scalar>& omega = mesh.lmodel()->omega();

    const scalar cs2 = mesh.lmodel()->cs2();

    const scalar q = mesh.lmodel()->q();


    
    // Initialize laplacian

    scalar lap(0);



    // Move over velocities
    
    for( uint k = 1 ; k < q ; k++ ) {

	int nbid = nb[id][k];

	if( nbid != -1 ) {

	    inverse  ?  lap += 2.0 * omega[k] * ((1/field[nbid] - 1/field[id]) ) / cs2  :  lap += 2.0 * omega[k] * (field[nbid] - field[id]) / cs2;	    

	}

	
	// If neighbour does not exist, use virtual node
	
	else {

	    
	    int first( mesh.vnode(id, reverse[k]) ),
		second( mesh.vnode(id, reverse[k], false) );

	    scalar fval(0);

	    if( second != -1 ) {

		inverse ? fval = 2/field[first] - 1/field[second] : fval = 2*field[first] - field[second];

	    }

	    else {

		inverse ? fval = 1/field[first] : fval = field[first];

	    }


	    // Add virtual node contribution

	    inverse  ?  lap += 2.0 * omega[k] * ((fval - 1/field[id]) ) / cs2  :  lap += 2.0 * omega[k] * (fval - field[id]) / cs2;	    	    


	    

	    // int otherNb = nb[id][reverse[k]];

	    // if( otherNb != -1 ) {

	    // 	scalar extpNb = 2.0*field[id]-field[otherNb];

	    // 	// inverse  ?  lap += 2.0 * omega[k]  * ((1.0/extpNb - 1/field[id]) ) / cs2  :  lap += 2.0 * omega[k] * (extpNb - field[id]) / cs2;

	    // 	if( inverse ) {

	    // 	    lap += 2.0 * omega[k]  * ((1.0/extpNb - 1/field[id]) ) / cs2;

	    // 	}

	    // 	else {

	    // 	    lap += 2.0 * omega[k] * (extpNb - field[id]) / cs2;

	    // 	}

	    // }

	    // else {



	    // }
	    

	}

	

    }


    return lap;
    

}









/** Field average */

const scalar scalarField::average() const {

    
    scalar avg(0);

   
    // Sum local nodes locally

    scalar localSum = 0;

    const uint np = mesh.local();

    for( uint i = 0 ; i < np ; i++ )
	localSum += field[i];

    


    // Apply collective reduction

    scalar globalSum = 0;

    int nelem = 0;

    int nlocal = mesh.local();
    
    if( mesh.wsize() > 1 ) {

	#ifdef DP
	    
	MPI_Allreduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	#elif SP

	MPI_Allreduce(&localSum, &globalSum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	#endif

	    

	MPI_Allreduce(&nlocal, &nelem, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);	
	    
    }

    else {
	    
	globalSum = localSum;

	nelem = mesh.local();

    }



    avg = globalSum / nelem;



    return avg;


}








/** Update field values from time index (not time value) */

const void scalarField::update( const uint& t ) {



    // MPI file pointer

    MPI_File file;

    MPI_File_open( MPI_COMM_WORLD,
		   ( "lattice." + name + "_" + to_string(t) ).c_str(),
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
