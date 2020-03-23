#include <iostream>

#include <fieldConstructor.H>

#include <initialScalarField.H>


using namespace std;



int main( int argc, char **argv ) {

    
    // Initialize mpi

    int pid, world;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD,&pid);

    MPI_Comm_size(MPI_COMM_WORLD,&world);



    if(pid == 0) {
    	cout << "                    " << endl;
    	cout << "     o-----o-----o  " << endl;
    	cout << "     | -   |   - |  " << endl;
    	cout << "     |   - | -   |  " << endl;
    	cout << "     o<----o---->o   Initial fields distribution" << endl;
    	cout << "     |   - | -   |  " << endl;
    	cout << "     | -   |   - |  " << endl;
    	cout << "     o-----o-----o  " << endl << endl;
    }




    // Lattice mesh creation
    
    latticeMesh mesh(pid);


    // Simulation handler

    timeOptions Time(pid);


    // Option dictionary

    dictionary dict("start/initialFields");

    
    // Set scalarFields

    vector<string> scfields = dict.bracedEntry("scalarFields");    

    for(auto fname : scfields) {

	if(pid == 0)
	    cout << "Initial values for scalar field " << fname << endl << endl; 
	

	// Field construction without reading

	scalarField field( mesh, Time, fname, IO::NO_READ, IO::MUST_WRITE );

	

	// Hand-coded, but...
	// Read types as list. Apply one ofter another

	vector<string> functionType = dict.bracedEntriesNames( fname + "/internalField" );

	initialScalarField itfield;

	for(auto ft : functionType)
	    itfield.updateField(field, mesh, fname, ft);

	

	
	// Write field

	field.write();

    }


    // Set vectorFields

    vector<string> vfields = dict.bracedEntry("vectorFields");

    for(auto fname : vfields) {


	if(pid == 0)
	    cout << "Initial values for vector field " << fname << endl << endl; 
	

	
	// Field construction without reading

	vectorField field( mesh, Time, fname, IO::NO_READ, IO::MUST_WRITE );


	// Write field

	field.write();

    }


    // Set pdfFields

    vector<string> pfields = dict.bracedEntry("pdfFields");

    for(auto fname : pfields) {


	if(pid == 0)
	    cout << "Initial values for pdf field " << fname << endl << endl; 
	

	// Field construction without reading

	pdfField field( mesh, Time, fname, IO::NO_READ, IO::MUST_WRITE );


	// Write field

	field.write();

    }    



    // Update case file

    Time.updateCaseFile();
    	
    

    MPI_Finalize();

}
