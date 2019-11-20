// #include <updateCaseFile.h>
// #include <vtkInfo.h>
// #include <readVTKInfo.h>
// #include <stdlib.h>
// #include <string.h>
// #include <stdio.h>

#include "updateCaseFile.H"

#include <cstring>

#include <iostream>

void updateCaseFile( latticeMesh_C& mesh ) {


    if(mesh.parallel.pid == 0) {


    
	FILE *outFile;

	uint ns = 0;

	uint *tsteps = NULL;
    


	// Try to open .case file to get time steps

	outFile = fopen("lbm.case", "r");


	if( outFile ) {


	    // Read Number of steps

	    char word[256];

	    uint find = 0;

	    uint status;
	    
	    while(   ( fscanf(outFile, "%s", word) != EOF )   &&   ( find == 0 )   ) {

		if( strcmp(word, "number") == 0 ) {

		    status = fscanf(outFile, "%s", word);

		    if(   ( status != EOF )   &&   ( strcmp(word, "of") == 0 ) ) {

			status = fscanf(outFile, "%s", word);

			if(   ( status != EOF )   &&   ( strcmp(word, "steps:") == 0 ) ) {
			    
			    status = fscanf(outFile, "%s", word);

			    ns = atoi( word );

			    tsteps = (uint*)malloc( (ns+1) * sizeof(uint) );

			    find = 1;
			
			}			
			
		    }

		}

	    }




	    // Read Time list

	    rewind( outFile );
	    
	    while(   ( fscanf(outFile, "%s", word) != EOF )   ) {

		if( strcmp(word, "time") == 0 ) {

		    status = fscanf(outFile, "%s", word);

		    if(   ( status != EOF )   &&   ( strcmp(word, "values:") == 0 ) ) {

			uint k;

			for( k = 0 ; k < ns ; k++ ) {

			    status = fscanf(outFile, "%s", word);

			    tsteps[k] = atoi( word );

			}
			
		    }

		}

	    }



	    // Add current time

	    // if( mesh.time.current != tsteps[ns-1] ) {

	    ns++;
		
	    tsteps[ns-1] = 0;

	    // }


	    fclose(outFile);
	    
	    

	}

	else {

	    ns = 1;

	    tsteps = (uint*)malloc( sizeof(uint) );

	    tsteps[0] = 0;

	}    







	// // .case file

	// vtkInfo vtk = readVTKInfo();

	

	// // Open File
	    
	// outFile = fopen("lbm.case", "w");


	// fprintf(outFile,"#\n");
	// fprintf(outFile,"# EnSight 7.4.1 ((n))\n");
	// fprintf(outFile,"# Case File: lattice.case\n");
	// fprintf(outFile,"\n");
	// fprintf(outFile,"FORMAT\n");
	// fprintf(outFile,"\n");
	// fprintf(outFile,"type:  ensight gold\n");
	// fprintf(outFile,"\n");
	// fprintf(outFile,"GEOMETRY\n");
	// fprintf(outFile,"\n");
	// fprintf(outFile,"model:                     lattice.geo\n");
	// fprintf(outFile,"\n");
	// fprintf(outFile,"VARIABLE\n");
	// fprintf(outFile,"\n");

	// uint fid;
		
	// for( fid = 0 ; fid < vtk.nscalar ; fid++ ) {

	//     fprintf(outFile,"scalar per node:           %s lattice.%s_*\n",vtk.scalarFields[fid],vtk.scalarFields[fid]);
	    
	// }

	// for( fid = 0 ; fid < vtk.npdf ; fid++ ) {

	//     uint k;

	//     for( k = 0 ; k < mesh.lattice.Q ; k++ ) {
	    
	// 	fprintf(outFile,"scalar per node:           %s%d lattice.%s%d_*\n",vtk.pdfFields[fid],k,vtk.pdfFields[fid],k);

	//     }
	    
	// }
	
	// for( fid = 0 ; fid < vtk.nvector ; fid++ ) {

	//     fprintf(outFile,"vector per node:           %s lattice.%s_*\n",vtk.vectorFields[fid],vtk.vectorFields[fid]);
	    
	// }


	
	// fprintf(outFile,"\n");
	
	// fprintf(outFile,"TIME\n");

	// fprintf(outFile,"time set:                  1\n");

	// fprintf(outFile,"number of steps:           %d\n", ns);

	// fprintf(outFile,"filename start number:     0\n");

	// fprintf(outFile,"filename increment:        1\n");

	// fprintf(outFile,"time values:               %d\n", tsteps[0]);

	// for( fid = 1 ; fid < ns ; fid++ ) {

	//     fprintf(outFile,"                           %d\n", tsteps[fid]);

	// }

	// fclose(outFile);



	// free(tsteps);


    }
    
}
