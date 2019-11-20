#include "updateCaseFile.H"

#include <iostream>

#include <dictionary.H>

#include <fstream>

using namespace std;


void updateCaseFile( latticeMesh_C& mesh ) {


    // Dictionary

    dictionary dict("start/initialFields");
    
    
    // Fields to be updated

    vector<string> scalarFields = dict.bracedEntry( "scalarFields" );

    vector<string> vectorFields = dict.bracedEntry( "vectorFields" );

    vector<string> pdfFields = dict.bracedEntry( "pdfFields" );

	

    // Open File

    ofstream outFile;

    outFile.open( "lbm.case" );


    // Header

    outFile << "#" << endl;
    outFile << "# EnSight 7.4.1 ((n))\n";
    outFile << "# Case File: lattice.case\n";
    outFile << "\n";
    outFile << "FORMAT\n";
    outFile << "\n";
    outFile << "type:  ensight gold\n";
    outFile << "\n";
    outFile << "GEOMETRY\n";
    outFile << "\n";
    outFile << "model:                     lattice.geo\n";
    outFile << "\n";
    outFile << "VARIABLE\n";
    outFile << "\n";


    // Scalar fields

    for( auto f : scalarFields )
    	outFile << "scalar per node:           " << f << " lattice." << f << "_*" << endl;


    // PDF fields

    for( auto f : pdfFields ) {

    	for( uint k = 0 ; k < mesh.mesh.Q ; k++ )	
    	    outFile << "scalar per node:           " << f << k << " lattice." << f << k << "_*" << endl;

    }


    // Vector fields

    for( auto f : vectorFields )
    	outFile << "vector per node:           " << f << " lattice." << f << "_*" << endl;
    


	
    outFile << "\n";	
    outFile << "TIME\n";
    outFile << "time set:                  1\n";
    outFile << "number of steps:           1\n";
    outFile << "filename start number:     0\n";
    outFile << "filename increment:        1\n";
    outFile << "time values:               0\n";



    outFile.close();

    
}
