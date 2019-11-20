#include "computeVirtualNodes.H"

using namespace std;

void computeVirtualNodes( const basicMesh& mesh, std::vector< std::vector<int> >& virtualNodes, latticeModel* lbmodel ) {


    // Lattice constants

    const vector< vector<int> > vel = lbmodel->lvel();

    vector<uint> reverse = lbmodel->reverse();

    const uint q = lbmodel->q();

    

    // Use boundary nodes only 
    
    vector<uint> dist( q );

    vector<uint> vid( q );



    // Compute maximum number of virtual nodes

    uint maxvirt = 0;

    for( uint bid = 0 ; bid < mesh.bd.nbd ; bid++ )
	maxvirt += mesh.bd.nbdelem[bid];

    maxvirt = maxvirt * mesh.Q;



    // First attempt
    
    vector< vector<int> > vnodes;

    vnodes.resize( maxvirt );

    for( uint i = 0 ; i < maxvirt ; i++ )
	vnodes[i].resize(4,-1);

	


    // Count total number of virtual

    uint nvirtual(0);


    
    

    // Move over boundary nodes 

    for( uint bid = 0 ; bid < mesh.bd.nbd ; bid++ ) {


	// Move over boundary elements
	    
	for( uint id = 0 ; id < mesh.bd.nbdelem[bid] ; id++ ) {


	    // Node id

	    uint node = mesh.bd.bdPoints[bid][id];

		
	    // Move over velocities and check for non-existing neighbour

	    for( uint k = 0 ; k < q ; k++ ) {
		    
		if(mesh.nb[node][ reverse[k] ] == -1) {


		    // First compute distance to other nodes conecting "node"

		    for( uint j = 0 ; j < mesh.Q ; j++ ) {

			dist[j] = (vel[j][0] - vel[k][0]) * (vel[j][0] - vel[k][0])
			        + (vel[j][1] - vel[k][1]) * (vel[j][1] - vel[k][1])
			        + (vel[j][2] - vel[k][2]) * (vel[j][2] - vel[k][2]);

			vid[j] = j;

		    }


		    // Sort distances

		    uint perm = 1;

		    while( perm != 0 ) {

		    	perm = 0;
			    
		    	for( uint j = 1 ; j < q ; j++ ) {

		    	    if( dist[j] < dist[j-1] ) {

		    		uint swpDist = dist[j-1],
		    		    swpVid = vid[j-1];

		    		dist[j-1] = dist[j];

		    		dist[j] = swpDist;


		    		vid[j-1] = vid[j];

		    		vid[j] = swpVid;

		    		perm++;

		    	    }

		    	}

		    }


			
		    // First distance is node. Check from second

		    uint d = dist[1];

		    uint maxv = 1;

		    for( uint j = 2 ; j < mesh.Q ; j++ ) {

		    	if( dist[j] == d )
		    	    maxv = j;

		    }



		    // Check if node exists at related velocity and asign closest wall node

		    int stj = 1, endj = maxv;

		    int wallNode = -1;

		    int wallVel = 0;

		    while(wallNode == -1) {

		    	for( int j = stj ; j <= endj ; j++ ) {

		    	    int neigh = mesh.nb[node][ reverse[vid[j]] ];
			    
		    	    if( neigh != -1 ) {

		    		wallNode = neigh;

		    		wallVel = vid[j];

		    	    }

		    	}


		    	// Move to next distance
			    
		    	if( wallNode == -1 ) {

		    	    stj = endj + 1;

		    	    d = dist[stj];

		    	    for( int j = stj ; j < (int)q ; j++ ) {

		    		if( dist[j] == d )
		    		    endj = j;

		    	    }
				

		    	}

		    }





		    // Finally detect fluid node in same direction

		    int sep[3] = { vel[wallVel][0] - vel[k][0], vel[wallVel][1] - vel[k][1], vel[wallVel][2] - vel[k][2] };

		    uint fvel = 0;
			
		    for( uint j = 0 ; j < q ; j++ ) {

		    	int sep2[3] = { vel[j][0] - vel[wallVel][0], vel[j][1] - vel[wallVel][1], vel[j][2] - vel[wallVel][2] };

		    	if(  ( sep2[0] == sep[0] )   &&   ( sep2[1] == sep[1] )   &&   ( sep2[2] == sep[2] )   )
		    	    fvel = j;

		    }


		    // Add to nodes

			
		    vnodes[nvirtual][0] = node;
		    vnodes[nvirtual][1] = k;
		    vnodes[nvirtual][2] = wallNode;
		    vnodes[nvirtual][3] = mesh.nb[node][reverse[fvel]];
		    nvirtual++;

		       

			
			

			

		}

	    }
		

	}
	 
	   	    

    }












    // Copy to reduced array
    
    virtualNodes.resize( nvirtual );

    for(uint i = 0 ; i < nvirtual ; i++)
	virtualNodes[i].resize(4, -1);
	

    for( uint id = 0 ; id < nvirtual ; id++ ) {

	for( uint k = 0 ; k < 4 ; k++ ) {

	    virtualNodes[id][k] = vnodes[id][k];

	}

    }
    


    
    
}





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
