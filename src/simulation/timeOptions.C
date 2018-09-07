#include <timeOptions.H>

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

    }
	
    
    return wrt;

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

	    uint ns;


	    
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

    	}    


	

	# warning Incomplete case file writing


    // 	// .case file

    // 	vtkInfo vtk = readVTKInfo();

	

    // 	// Open File
	    
    // 	cfile = fopen("lbm.case", "w");


    // 	fprintf(cfile,"#\n");
    // 	fprintf(cfile,"# EnSight 7.4.1 ((n))\n");
    // 	fprintf(cfile,"# Case File: lattice.case\n");
    // 	fprintf(cfile,"\n");
    // 	fprintf(cfile,"FORMAT\n");
    // 	fprintf(cfile,"\n");
    // 	fprintf(cfile,"type:  ensight gold\n");
    // 	fprintf(cfile,"\n");
    // 	fprintf(cfile,"GEOMETRY\n");
    // 	fprintf(cfile,"\n");
    // 	fprintf(cfile,"model:                     lattice.geo\n");
    // 	fprintf(cfile,"\n");
    // 	fprintf(cfile,"VARIABLE\n");
    // 	fprintf(cfile,"\n");

    // 	uint fid;
		
    // 	for( fid = 0 ; fid < vtk.nscalar ; fid++ ) {

    // 	    fprintf(cfile,"scalar per node:           %s lattice.%s_*\n",vtk.scalarFields[fid],vtk.scalarFields[fid]);
	    
    // 	}

    // 	for( fid = 0 ; fid < vtk.npdf ; fid++ ) {

    // 	    uint k;

    // 	    for( k = 0 ; k < mesh->lattice.Q ; k++ ) {
	    
    // 		fprintf(cfile,"scalar per node:           %s%d lattice.%s%d_*\n",vtk.pdfFields[fid],k,vtk.pdfFields[fid],k);

    // 	    }
	    
    // 	}
	
    // 	for( fid = 0 ; fid < vtk.nvector ; fid++ ) {

    // 	    fprintf(cfile,"vector per node:           %s lattice.%s_*\n",vtk.vectorFields[fid],vtk.vectorFields[fid]);
	    
    // 	}


	
    // 	fprintf(cfile,"\n");
	
    // 	fprintf(cfile,"TIME\n");

    // 	fprintf(cfile,"time set:                  1\n");

    // 	fprintf(cfile,"number of steps:           %d\n", ns);

    // 	fprintf(cfile,"filename start number:     0\n");

    // 	fprintf(cfile,"filename increment:        1\n");

    // 	fprintf(cfile,"time values:               %d\n", tsteps[0]);

    // 	for( fid = 1 ; fid < ns ; fid++ ) {

    // 	    fprintf(cfile,"                           %d\n", tsteps[fid]);

    // 	}

    // 	fclose(cfile);



    // 	free(tsteps);


    }


}
