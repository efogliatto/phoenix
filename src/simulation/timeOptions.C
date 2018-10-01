#include <timeOptions.H>

#include <mpi.h>

using namespace std;


/** Default constructor */

timeOptions::timeOptions(uint id) : pid(id), beg(clock_::now()) {


    dictionary dict("properties/simulation");


    
    start = (uint)dict.lookUp<int>("startTime");

    current = start;

    stp = 0;
    

    
    end = (uint)dict.lookUp<int>("endTime");

    writeInterval = (uint)dict.lookUp<int>("writeInterval");


    dataFormat = "binary";
  

}




/** Default destructor */

timeOptions::~timeOptions() {}




/** Seconds since class instantiation*/

const scalar timeOptions::elapsed() const {

    return (scalar)chrono::duration_cast<second_> ( clock_::now() - beg ).count();   

}



/** Update time: move to next time step */

const bool timeOptions::update() {

    
    bool upd = true;

    
    // Update time
    
    current++;
    
    stp++;


    if ( current > end )	
	upd = false;
	
    
    if( stp > writeInterval )	
	stp = 1;

	

    return upd;

}





/** Write flag */

const bool timeOptions::write() const {

    bool wrt(false);

    
    if ( stp == writeInterval )  {
	
	wrt = true;

	updateCaseFile();

	MPI_Barrier(MPI_COMM_WORLD);
	
    }
	
    
    return wrt;

}





/** Add simulation file names */

const void timeOptions::addScalarField( const string& name ) {

    scalarFields.push_back(name);
    
}

const void timeOptions::addVectorField( const string& name ) {

    vectorFields.push_back(name);

}

const void timeOptions::addPdfField( const string& name, const uint& q ) {

    pdfFields.push_back( make_tuple(name,q) );

}





/** Update ensight case file */

const void timeOptions::updateCaseFile() const {
    

    if( pid == 0 ) {
  

	// Time steps list
	
    	vector<uint> tsteps;    	
	

    	// Try to open .case file to get time steps

	ifstream cfile( "lbm.case" );

	if( cfile.is_open() ) {	


    	    // Read Number of steps

    	    string word;

    	    bool find(false);

	    uint ns(0);


	    
    	    while(  ( !cfile.eof() )   &&   ( find == false )   ) {

		cfile >> word;
		
    		if( word == "number" ) {

		    cfile >> word;

		    if (  ( !cfile.eof() )   &&   ( word == "of" )   ) {		    

			cfile >> word;

			if (  ( !cfile.eof() )   &&   ( word == "steps:" )   ) {		       
			    
			    cfile >> word;

    			    ns = stoi( word );

    			    find = true;
			
    			}
			
    		    }

    		}

    	    }




    	    // Read Time list

	    cfile.clear();
	    
	    cfile.seekg(0);
	    
    	    while( !cfile.eof() ) {

		cfile >> word;

    		if( word =="time" ) {

		    cfile >> word;		   

		    if(  (!cfile.eof())  &&  (word == "values:")  ) {

    			for( uint k = 0 ; k < ns ; k++ ) {

			    cfile >> word;			    

    			    tsteps.push_back( stoi( word ) );

    			}
			
    		    }

    		}

    	    }



    	    // Add current time

    	    if( current != tsteps.back() ) {
		
    		tsteps.push_back( current );

    	    }


    	    cfile.close();	    	    

    	}

    	else {

    	    tsteps.push_back(0);

    	    // cfile.close();
	    

    	}    

	

	

    	// Open File

	ofstream outFile("lbm.case");

	if( outFile.is_open() ) {

	    outFile << "#\n";
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
		
	    for( auto s : scalarFields )
		outFile << "scalar per node:           " << s << " lattice." << s << "_*" << endl;
	    

	
	    for( auto s : pdfFields ) {       

		for( uint k = 0 ; k < get<1>(s) ; k++ ) {

		    outFile << "scalar per node:           " << get<0>(s) << k << " lattice." << get<0>(s) << k << "_*" << endl;

		}
	    
	    }
	

	    for( auto s : vectorFields )
		outFile << "vector per node:           " << s << " lattice." << s << "_*" << endl;
	    

	
	    outFile << endl;
	
	    outFile << "TIME" << endl;

	    outFile << "time set:                  1" << endl;

	    outFile << "number of steps:           " << tsteps.size() << endl;

	    outFile << "filename start number:     0" << endl;

	    outFile << "filename increment:        1" << endl;

	    outFile << "time values:               " << tsteps[0] << endl;

	    for( uint fid = 1 ; fid < tsteps.size() ; fid++ )
		outFile << "                           " << tsteps[fid] << endl;



	
	    outFile.close();


	}


	else {

	    cout << " [ERROR]  Unable to open lbm.case " << endl << endl;

	    exit(1);

	}


    }


}









/** Match time to ensight index */

const uint timeOptions::timeToIndex( const uint& tid ) const {


    uint idx(0);

		

    // Try to open .case file to get time steps

    ifstream cfile( "lbm.case" );

    if( cfile.is_open() ) {	


	// Read Number of steps

	string word;

	bool find(false);

	uint ns(0);


	    
	while(  ( !cfile.eof() )   &&   ( find == false )   ) {

	    cfile >> word;
		
	    if( word == "number" ) {

		cfile >> word;

		if (  ( !cfile.eof() )   &&   ( word == "of" )   ) {		    

		    cfile >> word;

		    if (  ( !cfile.eof() )   &&   ( word == "steps:" )   ) {		       
			    
			cfile >> word;

			ns = stoi( word );

			find = true;
			
		    }
			
		}

	    }

	}




	// Read Time list

	cfile.clear();
	    
	cfile.seekg(0);
	    
	while( !cfile.eof() ) {

	    cfile >> word;

	    if( word =="time" ) {

		cfile >> word;		   

		if(  (!cfile.eof())  &&  (word == "values:")  ) {

		    for( uint k = 0 ; k < ns ; k++ ) {

			cfile >> word;

			if( stoi(word) == (int)tid )
			    idx = k;

		    }
			
		}

	    }

	}





	cfile.close();	    	    

    }
       

    

    return idx;
    

}
