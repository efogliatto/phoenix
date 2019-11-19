/*

  latticeMeshPartition

  Mesh subdivision for parallel processing

 */


#include <iostream>

#include "readBasicMesh.H"

#include <dictionary.H>

#include "kmetisDecomp.H"



/* #include <io.h> */
/* #include <dictIO.h> */
/* #include <latticeModel.h> */
/* #include <basic.h> */
/* #include <writeLatticeMesh.h> */


/* // Standard decomposition */
/* void standardDecomp( uint* owner, basicMesh* mesh, uint np ); */

/* // kmetis decomposition */
/* void kmetisDecomp( uint* owner, basicMesh* mesh, uint np ); */

/* // Local Indexing */
/* void localIndexing ( basicMesh* mesh, int** local, int** nGhosts, uint* owner, uint np ); */


using namespace std;


int main(int argc, char** argv) {



    cout << "                    " << endl;
    cout << "     o-----o-----o  " << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     |   - | -   |  latticeMeshPartition" << endl;
    cout << "     o<----o---->o  " << endl;
    cout << "     |   - | -   |   Mesh decomposition" << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     o-----o-----o  " << endl << endl;


    uint status;


    // Read full mesh

    basicMesh mesh = readBasicMesh();
    

    // Total number of processes

    dictionary parallelDict("properties/parallel");

    uint np( (uint)parallelDict.lookUpOrDefault<scalar>("numProc",1) );

   

    // Decomposition method
    
    string method( parallelDict.lookUpOrDefault<string>("method","standard") );
    



    /* int** virtualNodes; */

    /* uint nvirtual = 0; */



    /* // ******************************************************************** // */
    /* //                           Virtual nodes                              // */
    /* // ******************************************************************** // */

    /* { */


    /* 	printf("Virtual nodes\n\n"); */
	

    /* 	// Lattice model */

    /* 	latticeInfo lattice = setLatticeInfo(); */
	

	
    /* 	// Use boundary nodes only */
	
    /* 	uint bid, id, node, k; */

    /* 	uint* dist = (uint*)malloc( (mesh.Q) * sizeof(uint) ); */

    /* 	uint* vid = (uint*)malloc( (mesh.Q) * sizeof(vid) ); */



    /* 	// Compute maximun number of virtual nodes */

    /* 	uint maxvirt = 0; */

    /* 	for( bid = 0 ; bid < mesh.bd.nbd ; bid++ ) { */

    /* 	    maxvirt += mesh.bd.nbdelem[bid]; */

    /* 	} */

    /* 	maxvirt = maxvirt * mesh.Q; */

    /* 	int** vnodes = matrixIntAlloc( maxvirt, 4, -1); */


	
	
	

    /* 	// Move over boundary nodes */

    /* 	for( bid = 0 ; bid < mesh.bd.nbd ; bid++ ) { */


    /* 	    // Move over boundary elements */
	    
    /* 	    for( id = 0 ; id < mesh.bd.nbdelem[bid] ; id++ ) { */


    /* 		// Node id */

    /* 		node = mesh.bd.bdPoints[bid][id]; */

		
    /* 		// Move over velocities and check for non-existing neighbour */

    /* 		for( k = 0 ; k < mesh.Q ; k++ ) {		     */
		    
    /* 		    if(mesh.nb[node][ lattice.reverse[k] ] == -1) { */


    /* 			// First compute distance to other nodes conecting "node" */

    /* 			uint j; */

    /* 			for( j = 0 ; j < mesh.Q ; j++ ) { */

    /* 			    dist[j] = (lattice.vel[j][0] - lattice.vel[k][0]) * (lattice.vel[j][0] - lattice.vel[k][0]) */
    /* 				    + (lattice.vel[j][1] - lattice.vel[k][1]) * (lattice.vel[j][1] - lattice.vel[k][1]) */
    /* 				    + (lattice.vel[j][2] - lattice.vel[k][2]) * (lattice.vel[j][2] - lattice.vel[k][2]); */

    /* 			    vid[j] = j; */

    /* 			} */


    /* 			// Sort distances */

    /* 			uint perm = 1; */

    /* 			while( perm != 0 ) { */

    /* 			    perm = 0; */
			    
    /* 			    for( j = 1 ; j < mesh.Q ; j++ ) { */

    /* 				if( dist[j] < dist[j-1] ) { */

    /* 				    uint swpDist = dist[j-1], */
    /* 					swpVid = vid[j-1]; */

    /* 				    dist[j-1] = dist[j]; */

    /* 				    dist[j] = swpDist; */


    /* 				    vid[j-1] = vid[j]; */

    /* 				    vid[j] = swpVid; */

    /* 				    perm++; */

    /* 				} */

    /* 			    } */

    /* 			} */


			
    /* 			// First distance is node. Check from second */

    /* 			uint d = dist[1]; */

    /* 			uint maxv = 1; */

    /* 			for( j = 2 ; j < mesh.Q ; j++ ) { */

    /* 			    if( dist[j] == d ) */
    /* 				maxv = j; */

    /* 			} */



    /* 			// Check if node exists at related velocity and asign closest wall node */

    /* 			uint stj = 1, endj = maxv; */

    /* 			uint wallNode = -1; */

    /* 			uint wallVel = 0; */

    /* 			while(wallNode == -1) { */

    /* 			    for( j = stj ; j <= endj ; j++ ) { */

    /* 				uint neigh = mesh.nb[node][ lattice.reverse[vid[j]] ]; */
			    
    /* 				if( neigh != -1 ) { */

    /* 				    wallNode = neigh; */

    /* 				    wallVel = vid[j]; */

    /* 				} */

    /* 			    } */


    /* 			    // Move to next distance */
			    
    /* 			    if( wallNode == -1 ) { */

    /* 				stj = endj + 1; */

    /* 				d = dist[stj]; */

    /* 				for( j = stj ; j < mesh.Q ; j++ ) { */

    /* 				    if( dist[j] == d ) */
    /* 					endj = j; */

    /* 				}				 */
				

    /* 			    } */

    /* 			} */





    /* 			// Finally detect fluid node in same direction */

    /* 			int sep[3] = { lattice.vel[wallVel][0] - lattice.vel[k][0], lattice.vel[wallVel][1] - lattice.vel[k][1], lattice.vel[wallVel][2] - lattice.vel[k][2] }; */

    /* 			uint fvel = 0; */
			
    /* 			for( j = 0 ; j < mesh.Q ; j++ ) { */

    /* 			    int sep2[3] = { lattice.vel[j][0] - lattice.vel[wallVel][0], lattice.vel[j][1] - lattice.vel[wallVel][1], lattice.vel[j][2] - lattice.vel[wallVel][2] };			     */

    /* 			    if(  ( sep2[0] == sep[0] )   &&   ( sep2[1] == sep[1] )   &&   ( sep2[2] == sep[2] )   ) */
    /* 				fvel = j; */

    /* 			} */


    /* 			// Add to nodes */

    /* 			/\* printf("Node %d: %d %d %d\n", node, k, wallNode, mesh.nb[node][lattice.reverse[fvel]]); *\/ */
			
    /* 			vnodes[nvirtual][0] = node; */
    /* 			vnodes[nvirtual][1] = k; */
    /* 			vnodes[nvirtual][2] = wallNode; */
    /* 			vnodes[nvirtual][3] = mesh.nb[node][lattice.reverse[fvel]]; */
    /* 			nvirtual++; */

		       

			
			

			

    /* 		    } */

    /* 		} */
		

    /* 	    } */
	 
	   	    

    /* 	} */




    /* 	// Copy to reduced array */

    /* 	virtualNodes = matrixIntAlloc( nvirtual, 4, -1 ); */

    /* 	for( id = 0 ; id < nvirtual ; id++ ) { */

    /* 	    for( k = 0 ; k < 4 ; k++ ) { */

    /* 		virtualNodes[id][k] = vnodes[id][k]; */

    /* 	    } */

    /* 	} */
	    


	

    /* 	for( id = 0 ; id < maxvirt ; id++ ) */
    /* 	    free(vnodes[id]); */

    /* 	free(vnodes);	 */
	
    /* 	free(dist); */

    /* 	free(vid); */



    /* } */
    
    
    
    
    // ******************************************************************** //
    //                             Processors                               //
    // ******************************************************************** //


    cout << "Decomposing domain in " << np << " processors" << endl << endl;

    
    // Ownership array

    vector<uint> owner( mesh.nPoints );

    
    // Choose decomposition method
    
    if( method == "standard" ) {

    	// standardDecomp( owner, &mesh, np );

    }

    else {

    	if( method == "kmetis" ) {

    	    kmetisDecomp( owner, mesh, np );
	
    	}

    	else {

    	    cout << "\n\n  [ERROR]  Unable to recognize decomposition method " << method << "\n\n\n";

    	    exit(1);

    	}

    }





    /* // Resize local indices array */
    /* int** local   = matrixIntAlloc( mesh.nPoints, np, -1); */
    /* int** nGhosts = matrixIntAlloc( np, 2, -1); */
    /* /\* int** shared  = matrixIntAlloc( np, np, 0); *\/ */
    /* int** nrecv  = matrixIntAlloc( np, np, 0); */
    /* int** nsend  = matrixIntAlloc( np, np, 0); */

    
    /* // Creation of local indexing */
    /* localIndexing( &mesh, local, nGhosts, owner, np ); */


    /* // Total number of shared nodes */
    /* uint i, */
    /* 	 pid; */
    
    /* /\* for( i = 0 ; i < mesh.nPoints ; i++ ) { *\/ */

    /* /\* 	for( pid = 0 ; pid < np ; pid++) { *\/ */
	    
    /* /\* 	    if( local[i][pid] >= nGhosts[pid][0] ) { *\/ */

    /* /\* 		shared[pid][ owner[i] ]++; *\/ */

    /* /\* 	    } *\/ */

    /* /\* 	} *\/ */

    /* /\* } *\/ */



    /* // Individual send and receive ghosts */

    /* for( i = 0 ; i < mesh.nPoints ; i++ ) { */

    /* 	for( pid = 0 ; pid < np ; pid++) { */
	    
    /* 	    if( local[i][pid] >= nGhosts[pid][0] ) { */


    /* 		uint jj = owner[i]; */

    /* 		nrecv[jj][pid]++; */

    /* 		nsend[pid][jj]++; */
		

    /* 	    } */

    /* 	} */

    /* } */




    /* /\* for( i = 0 ; i < np ; i++ ) { *\/ */

    /* /\* 	for( pid = 0 ; pid < np ; pid++) { *\/ */

    /* /\* 	    printf("%d ", nrecv[i][pid]); *\/ */

    /* /\* 	} *\/ */

    /* /\* 	printf("\n"); *\/ */

    /* /\* } *\/ */
    
    /* /\* printf("\n\n"); *\/ */
    
    /* /\* for( i = 0 ; i < np ; i++ ) { *\/ */

    /* /\* 	for( pid = 0 ; pid < np ; pid++) { *\/ */

    /* /\* 	    printf("%d ", nsend[i][pid]); *\/ */

    /* /\* 	} *\/ */

    /* /\* 	printf("\n"); *\/ */

    /* /\* } *\/ */
    
    
    
    /* // Local mesh creation */
    /* latticeMesh* localMesh = (latticeMesh*)malloc( np * sizeof(latticeMesh) ); */

    /* // Move over meshes and look for recv ghosts */
    /* { */

	
    /* 	uint rpid,spid; */

    /* 	// Counter arrays */
    /* 	int** gcount = matrixIntAlloc( np, np, 0); */
	

    /* 	// Move over recv lattices. Basic info and resize arrays */
	
    /* 	for( rpid = 0 ; rpid < np ; rpid++ ) { */


    /* 	    // Add basic info */
	    
    /* 	    localMesh[rpid].parallel.pid = rpid; */

    /* 	    localMesh[rpid].parallel.worldSize = np; */

    /* 	    localMesh[rpid].parallel.nlocal = nGhosts[rpid][0]; */
	
    /* 	    localMesh[rpid].parallel.nghosts = nGhosts[rpid][1]; */




    /* 	    // Lattice model */

    /* 	    localMesh[rpid].lattice = setLatticeInfo(); */



    /* 	    /\* // Add sharing info and resize elements *\/ */
	    
    /* 	    /\* localMesh[rpid].parallel.shared = (uint*)malloc( np * sizeof(uint) ); *\/ */

    /* 	    /\* for( spid = 0 ; spid < np ; spid++ ) { *\/ */

    /* 	    /\* 	localMesh[rpid].parallel.shared[spid] = shared[rpid][spid]; *\/ */

    /* 	    /\* } *\/ */


    /* 	    // Add sharing info and resize elements */
	    
    /* 	    localMesh[rpid].parallel.nrg = (uint*)malloc( np * sizeof(uint) ); */

    /* 	    localMesh[rpid].parallel.nsg = (uint*)malloc( np * sizeof(uint) ); */

    /* 	    for( spid = 0 ; spid < np ; spid++ ) { */

    /* 		localMesh[rpid].parallel.nrg[spid] = nrecv[spid][rpid]; */

    /* 		localMesh[rpid].parallel.nsg[spid] = nsend[spid][rpid]; */

    /* 	    } */


    /* 	    // Resize ghost info */
	    
    /* 	    localMesh[rpid].parallel.recvGhosts = (uint**)malloc( np * sizeof(uint*) ); */

    /* 	    localMesh[rpid].parallel.sendGhosts = (uint**)malloc( np * sizeof(uint*) ); */
	    
    /* 	    for( spid = 0 ; spid < np ; spid++ ) { */

    /* 		localMesh[rpid].parallel.recvGhosts[spid] = (uint*)malloc( nrecv[spid][rpid] * sizeof(uint) ); */

    /* 		localMesh[rpid].parallel.sendGhosts[spid] = (uint*)malloc( nsend[spid][rpid] * sizeof(uint) ); */
		
    /* 	    } */




    /* 	    // Add points. First local, then ghost */

    /* 	    localMesh[rpid].mesh.nPoints = localMesh[rpid].parallel.nlocal + localMesh[rpid].parallel.nghosts; */
	    
    /* 	    localMesh[rpid].mesh.points = matrixIntAlloc( localMesh[rpid].mesh.nPoints, 3, 0 ); */

    /* 	    for( i = 0 ; i < mesh.nPoints ; i++ ) { */

    /* 		int lid = local[i][rpid]; */
		
    /* 		if( lid != -1 ) { */

    /* 		    localMesh[rpid].mesh.points[lid][0] = mesh.points[i][0]; */
    /* 		    localMesh[rpid].mesh.points[lid][1] = mesh.points[i][1]; */
    /* 		    localMesh[rpid].mesh.points[lid][2] = mesh.points[i][2]; */
		    
    /* 		} */

    /* 	    } */



    /* 	    // Add Neighbours. */
	    
    /* 	    localMesh[rpid].mesh.nb = matrixIntAlloc( localMesh[rpid].parallel.nlocal, mesh.Q, -1 ); */

    /* 	    localMesh[rpid].mesh.Q = mesh.Q; */

    /* 	    localMesh[rpid].lattice.Q = mesh.Q; */

    /* 	    for( i = 0 ; i < mesh.nPoints ; i++ ) { */

    /* 		int lid = local[i][rpid]; */
		
    /* 		if( lid < localMesh[rpid].parallel.nlocal ) { */


    /* 		    // Move over velocities */

    /* 		    uint velId; */
    /* 		    int nbid; */

    /* 		    for( velId = 0 ; velId < mesh.Q ; velId++ ) { */

    /* 			nbid = mesh.nb[i][velId]; */

    /* 			if( nbid != -1 ) { */
			    
    /* 			    localMesh[rpid].mesh.nb[lid][velId] = local[nbid][rpid]; */

    /* 			} */

    /* 		    } */
		    
		    
    /* 		} */

    /* 	    } */







    /* 	    // Add vtkCells */

    /* 	    localMesh[rpid].mesh.ncells = 0; */
	    
    /* 	    for( i = 0 ; i < mesh.ncells ; i++ ) { */

    /* 		uint cid, */
    /* 		     find = 0; */

    /* 		for( cid = 0 ; cid < mesh.cellType ; cid++ ) { */

    /* 		    // Check if all members are local */
    /* 		    if( local[ mesh.vtkCells[i][cid] ][rpid] == -1 ) { */

    /* 			find++; */

    /* 		    } */

    /* 		} */

    /* 		if( find == 0 ) { */

    /* 		    localMesh[rpid].mesh.ncells++; */

    /* 		} */

    /* 	    } */



    /* 	    // Resize and add */

    /* 	    localMesh[rpid].mesh.cellType = mesh.cellType; */
	    
    /* 	    localMesh[rpid].mesh.vtkCells = matrixIntAlloc( localMesh[rpid].mesh.ncells, mesh.cellType, -1); */

    /* 	    uint count = 0; */
	    
    /* 	    for( i = 0 ; i < mesh.ncells ; i++ ) { */

    /* 		uint cid, */
    /* 		     find = 0; */

    /* 		for( cid = 0 ; cid < mesh.cellType ; cid++ ) { */

    /* 		    // Check if all members are local */
    /* 		    if( local[ mesh.vtkCells[i][cid] ][rpid] == -1 ) { */

    /* 			find++; */

    /* 		    } */

    /* 		} */

    /* 		if( find == 0 ) { */

		    
    /* 		    for( cid = 0 ; cid < mesh.cellType ; cid++ ) { */

    /* 			localMesh[rpid].mesh.vtkCells[count][cid] = local[ mesh.vtkCells[i][cid] ][rpid]; */

    /* 		    } */
		    
		    
    /* 		    count++; */

    /* 		} */

    /* 	    } */








    /* 	    // Boundaries. Assign boundaries from original mesh */

    /* 	    localMesh[rpid].mesh.bd.nbd = mesh.bd.nbd; */
	    
    /* 	    localMesh[rpid].mesh.bd.nbdelem = (uint*)malloc( mesh.bd.nbd * sizeof(uint) ); */

    /* 	    localMesh[rpid].mesh.bd.bdPoints = (uint**)malloc( mesh.bd.nbd * sizeof(uint*) ); */

    /* 	    for( i = 0 ; i < localMesh[rpid].mesh.bd.nbd ; i++ ) { */

    /* 		localMesh[rpid].mesh.bd.nbdelem[i] = 0; */

    /* 		sprintf( localMesh[rpid].mesh.bd.bdNames[i], "%s", mesh.bd.bdNames[i] ); */

    /* 		uint bdpid; */

    /* 		for( bdpid = 0 ; bdpid < mesh.bd.nbdelem[i] ; bdpid++ ) { */

    /* 		    if( local[ mesh.bd.bdPoints[i][bdpid] ][rpid] < localMesh[rpid].parallel.nlocal ) { */

    /* 			localMesh[rpid].mesh.bd.nbdelem[i]++; */

    /* 		    } */

    /* 		} */

    /* 	    } */


	    
	    
	    
    /* 	    for( i = 0 ; i < localMesh[rpid].mesh.bd.nbd ; i++ ) { */

    /* 		count = 0; */
		
    /* 	    	localMesh[rpid].mesh.bd.bdPoints[i] = (uint*)malloc( localMesh[rpid].mesh.bd.nbdelem[i] * sizeof(uint) ); */

    /* 	    	uint bdpid; */

    /* 	    	for( bdpid = 0 ; bdpid < mesh.bd.nbdelem[i] ; bdpid++ ) { */

    /* 	    	    if( local[ mesh.bd.bdPoints[i][bdpid] ][rpid] < localMesh[rpid].parallel.nlocal ) { */

    /* 	    		localMesh[rpid].mesh.bd.bdPoints[i][count] = local[ mesh.bd.bdPoints[i][bdpid] ][rpid]; */

    /* 	    		count++; */

    /* 	    	    } */

    /* 	    	} */

    /* 	    } */




	    
	    
	    

    /* 	} */



	





    /* 	// Move over local lattices and add parallel info */

    /* 	for( i = 0 ; i < mesh.nPoints ; i++ ) { */
	
    /* 	    for( rpid = 0 ; rpid < np ; rpid++ ) { */

    /* 		if( local[i][rpid] >= nGhosts[rpid][0] ) { */
		    

    /* 		    // Add local index as recv ghost */

    /* 		    spid = owner[i]; */
		    
    /* 		    localMesh[rpid].parallel.recvGhosts[ spid ][ gcount[rpid][spid] ]   =  local[i][rpid]; */


    /* 		    // Add local index as send ghost */
		    
    /* 		    localMesh[spid].parallel.sendGhosts[ rpid ][ gcount[rpid][spid] ]   =  local[i][spid]; */

    /* 		    gcount[rpid][spid]++; */

    /* 		} */
	    
    /* 	    } */

    /* 	} */






	
	
    /* } */




    


    /* // Write lattice meshes */
    /* { */

    /* 	int status = system( "rm -rf processor*" ); */

    /* 	if (!status) { */

    /* 	    for( i = 0 ; i < np ; i++ ) { */

    /* 		writeLatticeMesh( &localMesh[i] ); */




    /* 		// Count number of virtual nodes per patch */

    /* 		uint count = 0; */

    /* 		uint id; */

    /* 		for( id = 0 ; id < nvirtual ; id++ ) { */

    /* 		    if( local[ virtualNodes[id][0] ][ i ] != -1 ) */
    /* 			count++; */

    /* 		} */
		

    /* 		// Write virtual nodes */

    /* 		char fname[100]; */

    /* 		FILE *outFile; */

    /* 		sprintf(fname,"processor%d/lattice/virtualNodes", i); */
    
    /* 		outFile = fopen(fname,"w"); */


		
    /* 		fprintf(outFile,"%d\n",count); */

    /* 		for( id = 0 ; id < nvirtual ; id++ ) { */

    /* 		    /\* if( local[ virtualNodes[id][0] ][ i ] != -1 ){ *\/ */
    /* 		    /\* 	fprintf(outFile,"%d %d %d %d\n", local[ virtualNodes[id][0] ][ i ], virtualNodes[id][1], local[ virtualNodes[id][2] ][ i ], local[ virtualNodes[id][3] ][ i ]); *\/ */
    /* 		    /\* } *\/ */

    /* 		} */

		
    /* 		fclose(outFile); */
		

    /* 	    } */

    /* 	} */
    /* } */

    

    



    /* printf("Finished domain decomposition\n\n"); */






    /* uint id; */
    
    /* for( id = 0 ; id < nvirtual ; id++ ) */
    /* 	free(virtualNodes[id]); */

    /* free(virtualNodes); */
    
    
    return 0;

}
