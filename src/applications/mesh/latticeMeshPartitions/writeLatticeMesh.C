#include "writeLatticeMesh.H"

#include "writeMeshToEnsight.H"

#include <iostream>

void writeLatticeMesh( latticeMesh_C& mesh ) {



    // Write basic lattice in Ensight format
    
    writeMeshToEnsight( mesh );


    
    uint i,j;

    // Write files
    char command[100];
    sprintf(command,"mkdir -p processor%d/lattice", mesh.parallel.pid);
    
    int status = system( command );
    if(status) {}


    FILE *outFile;

    
    // Local points

    sprintf(command,"processor%d/lattice/points", mesh.parallel.pid);
    
    outFile = fopen(command,"w");

    fprintf(outFile,"%d\n",mesh.parallel.nlocal);

    for( i = 0 ; i < mesh.parallel.nlocal ; i++ ) {

    	fprintf(outFile,"%d %d %d\n",mesh.mesh.points[i][0],mesh.mesh.points[i][1],mesh.mesh.points[i][2]);

    }
    
    fclose(outFile);



    // Ghost points

    sprintf(command,"processor%d/lattice/ghosts", mesh.parallel.pid);
    
    outFile = fopen(command,"w");

    fprintf(outFile,"%d\n",mesh.parallel.nghosts);

    for( i = mesh.parallel.nlocal ; i < mesh.mesh.nPoints ; i++ ) {

    	fprintf(outFile,"%d %d %d\n",mesh.mesh.points[i][0],mesh.mesh.points[i][1],mesh.mesh.points[i][2]);

    }
    
    fclose(outFile);


    
    
    // Neighbours
    
    sprintf(command,"processor%d/lattice/neighbours", mesh.parallel.pid);
    
    outFile = fopen(command,"w");

    uint dd;

    if(mesh.mesh.Q == 9) {dd = 2;}

    if(mesh.mesh.Q == 15) {dd = 3;}
    
    fprintf(outFile,"%d %d %d\n", mesh.parallel.nlocal, dd, mesh.mesh.Q);
    
    
    for( i = 0 ; i < mesh.parallel.nlocal ; i++ ) {

    	for( j = 0 ; j < mesh.mesh.Q ; j++ ) {

    	    fprintf(outFile,"%d ",mesh.mesh.nb[i][j]);

    	}

    	fprintf(outFile,"\n");

    }

    fclose(outFile);




    // Recv ghosts
    
    sprintf(command,"processor%d/lattice/recvGhosts", mesh.parallel.pid);
    
    outFile = fopen(command,"w");

    fprintf(outFile,"%d\n\n",mesh.parallel.worldSize);

    uint pid;

    for( pid = 0 ; pid < mesh.parallel.worldSize ; pid++ ) {

    	fprintf(outFile,"%d\n%d\n",pid,mesh.parallel.nrg[pid]);

    	for( i = 0 ; i < mesh.parallel.nrg[pid] ; i++ ) {

    	    fprintf(outFile,"%d\n",mesh.parallel.recvGhosts[pid][i]);

    	}

    	fprintf(outFile,"\n");

    }

    fclose(outFile);



    // Send ghosts
    
    sprintf(command,"processor%d/lattice/sendGhosts", mesh.parallel.pid);
    
    outFile = fopen(command,"w");

    fprintf(outFile,"%d\n\n",mesh.parallel.worldSize);

    for( pid = 0 ; pid < mesh.parallel.worldSize ; pid++ ) {

    	fprintf(outFile,"%d\n%d\n",pid,mesh.parallel.nsg[pid]);

    	for( i = 0 ; i < mesh.parallel.nsg[pid] ; i++ ) {

    	    fprintf(outFile,"%d\n",mesh.parallel.sendGhosts[pid][i]);

    	}

    	fprintf(outFile,"\n");

    }

    fclose(outFile);

    

    // Boundary

    sprintf(command,"processor%d/lattice/boundary", mesh.parallel.pid);
    
    outFile = fopen(command,"w");
	
    fprintf(outFile,"%d\n\n", mesh.mesh.bd.nbd);

    for( i = 0 ; i < mesh.mesh.bd.nbd ; i++ ) {

    	fprintf(outFile,"%s\n",mesh.mesh.bd.bdNames[i].c_str());
	
    	fprintf(outFile,"%d\n",mesh.mesh.bd.nbdelem[i]);

    	for( j = 0 ; j < mesh.mesh.bd.nbdelem[i] ; j++ ) {

    	    fprintf(outFile,"%d\n",mesh.mesh.bd.bdPoints[i][j]);

    	}

    	fprintf(outFile,"\n");
	

    }
	
	
    fclose(outFile);





    // Write VTK cells

    sprintf(command,"processor%d/lattice/vtkCells", mesh.parallel.pid);
    
    outFile = fopen(command,"w");
    
    fprintf(outFile,"%d %d\n", mesh.mesh.ncells, mesh.mesh.cellType);
    
    for( i = 0 ; i < mesh.mesh.ncells ; i++ ) {
	
    	for( j = 0 ; j < mesh.mesh.cellType ; j++ ) {
	    
    	    fprintf(outFile,"%d ",mesh.mesh.vtkCells[i][j]);
	    
    	}
	
    	fprintf(outFile,"\n");
	
    }
    
    fclose(outFile);
    

  


    
}
