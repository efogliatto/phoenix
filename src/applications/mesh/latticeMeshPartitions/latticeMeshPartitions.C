/*

  latticeMeshPartition

  Mesh subdivision for parallel processing

 */


#include <iostream>

#include "readBasicMesh.H"

#include <dictionary.H>

#include "kmetisDecomp.H"

#include "standardDecomp.H"

#include "localIndexing.H"

#include "../meshInclude/latticeMesh_C.H"

#include <latticeModelCreator.H>

#include "writeLatticeMesh.H"

#include "computeVirtualNodes.H"



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



    // Read full mesh

    basicMesh mesh = readBasicMesh();
    

    
    // Total number of processes

    dictionary parallelDict("properties/parallel");

    uint np( (uint)parallelDict.lookUpOrDefault<scalar>("numProc",1) );

   

    // Decomposition method
    
    string method( parallelDict.lookUpOrDefault<string>("method","standard") );


    
    // Read model name. Use D2Q9 as default

    dictionary ldict("properties/latticeProperties");
    
    latticeModelCreator lbm;
	
    latticeModel* lbmodel = lbm.create(ldict.lookUpOrDefault<string>("LBModel","D2Q9"));
    
    





    // ******************************************************************** //
    //                           Virtual nodes                              //
    // ******************************************************************** //

    cout << "Virtual nodes" << endl << endl;


    vector< vector<int> > virtualNodes;

    computeVirtualNodes( mesh, virtualNodes, lbmodel );

    // for(auto vn : virtualNodes) {

    // 	for( auto v : vn) {

    // 	    cout << v << " ";

    // 	}

    // 	cout << endl;

    // }
	
    
        
    
    
    
    // ******************************************************************** //
    //                             Processors                               //
    // ******************************************************************** //


    cout << "Decomposing domain in " << np << " processors" << endl << endl;

    
    // Ownership array

    vector<uint> owner( mesh.nPoints );

    
    // Choose decomposition method
    
    if( method == "standard" ) {

    	standardDecomp( owner, mesh, np );

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




    // Resize local indices array

    vector< vector<int> > local;

    local.resize(mesh.nPoints);

    for(uint i = 0 ; i < mesh.nPoints ; i++) {

	local[i].resize(np);

	std::fill( local[i].begin(),local[i].end(),-1 );

    }


    vector< vector<int> > nGhosts;

    nGhosts.resize(np);

    for(uint i = 0 ; i < np ; i++) {

	nGhosts[i].resize(2);

	std::fill( nGhosts[i].begin(),nGhosts[i].end(),-1 );

    }
    

    vector< vector<int> > nrecv;

    nrecv.resize(np);

    for(uint i = 0 ; i < np ; i++) {

	nrecv[i].resize(np);

	std::fill( nrecv[i].begin(),nrecv[i].end(), 0 );

    }


    vector< vector<int> > nsend;

    nsend.resize(np);

    for(uint i = 0 ; i < np ; i++) {

	nsend[i].resize(np);

	std::fill( nsend[i].begin(),nsend[i].end(), 0 );

    }
    
    


    
    // Creation of local indexing
    
    localIndexing( mesh, local, nGhosts, owner, np );
    



    // Individual send and receive ghosts

    for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	for( uint pid = 0 ; pid < np ; pid++) {
	    
    	    if( local[i][pid] >= nGhosts[pid][0] ) {


    		uint jj = owner[i];

    		nrecv[jj][pid]++;

    		nsend[pid][jj]++;
		

    	    }

    	}

    }




    
    
    
    // Local mesh creation

    vector<latticeMesh_C> localMesh( np );

    // Move over meshes and look for recv ghosts
    {

       

    	// Counter arrays

	vector< vector<int> > gcount;

	gcount.resize(np);

	for(uint i = 0 ; i < np ; i++) {

	    gcount[i].resize(np);

	    std::fill( gcount[i].begin(),gcount[i].end(), 0 );

	}
	
	
	

    	// Move over recv lattices. Basic info and resize arrays
	
    	for( uint rpid = 0 ; rpid < np ; rpid++ ) {


    	    // Add basic info
	    
    	    localMesh[rpid].parallel.pid = rpid;

    	    localMesh[rpid].parallel.worldSize = np;

    	    localMesh[rpid].parallel.nlocal = nGhosts[rpid][0];
	
    	    localMesh[rpid].parallel.nghosts = nGhosts[rpid][1];

	    localMesh[rpid].mesh.bd.bdNames.resize( mesh.bd.bdNames.size() );

	    localMesh[rpid].lattice_D = lbmodel->d();

	    localMesh[rpid].mesh.Q = lbmodel->q();
		
	    


    	    // // Lattice model

    	    // localMesh[rpid].lattice = setLatticeInfo();



    	    // Add sharing info and resize elements
	    
    	    localMesh[rpid].parallel.nrg.resize( np );

    	    localMesh[rpid].parallel.nsg.resize( np );

    	    for( uint spid = 0 ; spid < np ; spid++ ) {

    		localMesh[rpid].parallel.nrg[spid] = nrecv[spid][rpid];

    		localMesh[rpid].parallel.nsg[spid] = nsend[spid][rpid];

    	    }


    	    // Resize ghost info
	    
    	    localMesh[rpid].parallel.recvGhosts.resize( np );

    	    localMesh[rpid].parallel.sendGhosts.resize( np );
	    
    	    for( uint spid = 0 ; spid < np ; spid++ ) {

    		localMesh[rpid].parallel.recvGhosts[spid].resize( nrecv[spid][rpid] );

    		localMesh[rpid].parallel.sendGhosts[spid].resize( nsend[spid][rpid] );
		
    	    }




    	    // Add points. First local, then ghost

    	    localMesh[rpid].mesh.nPoints = localMesh[rpid].parallel.nlocal + localMesh[rpid].parallel.nghosts;
	    
    	    // localMesh[rpid].mesh.points = matrixIntAlloc( localMesh[rpid].mesh.nPoints, 3, 0 );
	    localMesh[rpid].mesh.points.resize( localMesh[rpid].mesh.nPoints );

	    for( uint i = 0 ; i < localMesh[rpid].mesh.points.size() ; i++ ) {

		localMesh[rpid].mesh.points[i].resize(3,0);

	    }
	    

    	    for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    		int lid = local[i][rpid];
		
    		if( lid != -1 ) {

    		    localMesh[rpid].mesh.points[lid][0] = mesh.points[i][0];
    		    localMesh[rpid].mesh.points[lid][1] = mesh.points[i][1];
    		    localMesh[rpid].mesh.points[lid][2] = mesh.points[i][2];
		    
    		}

    	    }



    	    // Add Neighbours.
	    
    	    // localMesh[rpid].mesh.nb = matrixIntAlloc( localMesh[rpid].parallel.nlocal, mesh.Q, -1 );
	    localMesh[rpid].mesh.nb.resize( localMesh[rpid].parallel.nlocal );

	    for( uint i = 0 ; i < localMesh[rpid].mesh.nb.size() ; i++ )
		localMesh[rpid].mesh.nb[i].resize(mesh.Q,-1);
	    

    	    localMesh[rpid].mesh.Q = mesh.Q;

    	    // // localMesh[rpid].lattice.Q = mesh.Q;

    	    for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	    	int lid = local[i][rpid];
		
    	    	if( lid < localMesh[rpid].parallel.nlocal ) {


    	    	    // Move over velocities

    	    	    for( uint velId = 0 ; velId < mesh.Q ; velId++ ) {

    	    		int nbid = mesh.nb[i][velId];

    	    		if( nbid != -1 ) {
			    
    	    		    localMesh[rpid].mesh.nb[lid][velId] = local[nbid][rpid];

    	    		}

    	    	    }
		    
		    
    	    	}

    	    }







    	    // Add vtkCells

    	    localMesh[rpid].mesh.ncells = 0;
	    
    	    for( uint i = 0 ; i < mesh.ncells ; i++ ) {

    		uint find = 0;

    		for( uint cid = 0 ; cid < mesh.cellType ; cid++ ) {

    		    // Check if all members are local
    		    if( local[ mesh.vtkCells[i][cid] ][rpid] == -1 ) {

    			find++;

    		    }

    		}

    		if( find == 0 ) {

    		    localMesh[rpid].mesh.ncells++;

    		}

    	    }



    	    // Resize and add

    	    localMesh[rpid].mesh.cellType = mesh.cellType;
	    
    	    // localMesh[rpid].mesh.vtkCells = matrixIntAlloc( localMesh[rpid].mesh.ncells, mesh.cellType, -1);
	    localMesh[rpid].mesh.vtkCells.resize( localMesh[rpid].mesh.ncells );

	    for( uint i = 0 ; i < localMesh[rpid].mesh.vtkCells.size() ; i++ )
		localMesh[rpid].mesh.vtkCells[i].resize(mesh.cellType, -1 );

    	    uint count = 0;
	    
    	    for( uint i = 0 ; i < mesh.ncells ; i++ ) {

    		uint find = 0;

    		for( uint cid = 0 ; cid < mesh.cellType ; cid++ ) {

    		    // Check if all members are local
    		    if( local[ mesh.vtkCells[i][cid] ][rpid] == -1 ) {

    			find++;

    		    }

    		}

    		if( find == 0 ) {

		    
    		    for( uint cid = 0 ; cid < mesh.cellType ; cid++ ) {

    			localMesh[rpid].mesh.vtkCells[count][cid] = local[ mesh.vtkCells[i][cid] ][rpid];

    		    }
		    
		    
    		    count++;

    		}

    	    }








    	    // Boundaries. Assign boundaries from original mesh

    	    localMesh[rpid].mesh.bd.nbd = mesh.bd.nbd;
	    
    	    localMesh[rpid].mesh.bd.nbdelem.resize( mesh.bd.nbd );

    	    localMesh[rpid].mesh.bd.bdPoints.resize( mesh.bd.nbd );

    	    for( uint i = 0 ; i < localMesh[rpid].mesh.bd.nbd ; i++ ) {

    	    	localMesh[rpid].mesh.bd.nbdelem[i] = 0;

    	    	localMesh[rpid].mesh.bd.bdNames[i] = mesh.bd.bdNames[i];

    	    	for( uint bdpid = 0 ; bdpid < mesh.bd.nbdelem[i] ; bdpid++ ) {

    	    	    if( local[ mesh.bd.bdPoints[i][bdpid] ][rpid] < localMesh[rpid].parallel.nlocal ) {

    	    		localMesh[rpid].mesh.bd.nbdelem[i]++;

    	    	    }

    	    	}

    	    }


	    
	    
	    
    	    for( uint i = 0 ; i < localMesh[rpid].mesh.bd.nbd ; i++ ) {

    		uint count = 0;
		
    	    	localMesh[rpid].mesh.bd.bdPoints[i].resize( localMesh[rpid].mesh.bd.nbdelem[i] );

    	    	for( uint bdpid = 0 ; bdpid < mesh.bd.nbdelem[i] ; bdpid++ ) {

    	    	    if( local[ mesh.bd.bdPoints[i][bdpid] ][rpid] < localMesh[rpid].parallel.nlocal ) {

    	    		localMesh[rpid].mesh.bd.bdPoints[i][count] = local[ mesh.bd.bdPoints[i][bdpid] ][rpid];

    	    		count++;

    	    	    }

    	    	}

    	    }




	    
	    
	    

    	}



	





    	// Move over local lattices and add parallel info

    	for( uint i = 0 ; i < mesh.nPoints ; i++ ) {
	
    	    for( uint rpid = 0 ; rpid < np ; rpid++ ) {

    		if( local[i][rpid] >= nGhosts[rpid][0] ) {
		    

    		    // Add local index as recv ghost

    		    uint spid = owner[i];
		    
    		    localMesh[rpid].parallel.recvGhosts[ spid ][ gcount[rpid][spid] ]   =  local[i][rpid];


    		    // Add local index as send ghost
		    
    		    localMesh[spid].parallel.sendGhosts[ rpid ][ gcount[rpid][spid] ]   =  local[i][spid];

    		    gcount[rpid][spid]++;

    		}
	    
    	    }

    	}






	
	
    }




    


    // Write lattice meshes
    {

    	int status = system( "rm -rf processor*" );

    	if (!status) {

    	    for( uint i = 0 ; i < np ; i++ ) {

    		writeLatticeMesh( localMesh[i] );




    		// Count number of virtual nodes per patch

    		uint count = 0;

    		for( uint id = 0 ; id < virtualNodes.size() ; id++ ) {
		    
		    if( local[ virtualNodes[id][0] ][ i ] != -1 )
			count++;

    		}
		

    		// Write virtual nodes

    		char fname[100];

    		FILE *outFile;

    		sprintf(fname,"processor%d/lattice/virtualNodes", i);
    
    		outFile = fopen(fname,"w");


		
    		fprintf(outFile,"%d\n",count);


		// for( uint id = 0 ; id < virtualNodes.size() ; id++ ) {

		//     if( local[ virtualNodes[id][0] ][ i ] != -1 ){
			
		// 	// fprintf(outFile,"%d %d %d %d\n", local[ virtualNodes[id][0] ][ i ], virtualNodes[id][1], local[ virtualNodes[id][2] ][ i ], local[ virtualNodes[id][3] ][ i ]);

		// 	fprintf(outFile,"%d %d ", local[ virtualNodes[id][0] ][ i ], virtualNodes[id][1]);

		// 	if( virtualNodes[id][2] != -1 ) {
			    
		// 	    fprintf(outFile,"%d ", local[ virtualNodes[id][2] ][ i ]);
			    
		// 	    if( virtualNodes[id][3] != -1 ) {
			    
		// 		fprintf(outFile,"%d\n", local[ virtualNodes[id][3] ][ i ]);			   
			    
		// 	    }
			    
		// 	}
			    
		//     }


		// }

		
    		fclose(outFile);
		

    	    }

    	}
    }

    

    



    cout << "Finished domain decomposition" << endl << endl;

   
    
    return 0;

}
