#include "computeVirtualNodes.H"

// #include <algorithm>

using namespace std;



// Auxiliary sort function

bool compareDist(vector<int>& d1, vector<int>& d2) {   return (d1[0] < d2[0]);  } 


void computeVirtualNodes( const basicMesh& mesh, std::vector< std::vector<int> >& vnodes, latticeModel* lbmodel ) {


    // Lattice constants

    const vector< vector<int> > vel = lbmodel->lvel();

    vector<uint> reverse = lbmodel->reverse();

    const uint q = lbmodel->q();

    

    // // Use boundary nodes only 
    // // distance[i] = distance, node
    
    // vector< vector<uint> > distance;

    // distance.resize(q*q);

    // for(uint i = 0 ; i < q ; i++)
    // 	distance[i].resize(2,-1);

    
    
    // vector<uint> nodeIds( q * q );

    // vector<uint> vid( q * q);



    // // Compute maximum number of virtual nodes

    // uint maxvirt = 0;

    // for( uint bid = 0 ; bid < mesh.bd.nbd ; bid++ )
    // 	maxvirt += mesh.bd.nbdelem[bid];

    // maxvirt = maxvirt * mesh.Q * mesh.Q;



    // First attempt
    
    // vector< vector<int> > vnodes;

    // vnodes.resize( maxvirt );

    // for( uint i = 0 ; i < maxvirt ; i++ )
    // 	vnodes[i].resize(4,-1);

	


    // // Count total number of virtual

    // uint nvirtual(0);


    
    

    // Move over boundary nodes 

    for( uint bid = 0 ; bid < mesh.bd.nbd ; bid++ ) {


	// Move over boundary elements
	    
	for( uint id = 0 ; id < mesh.bd.nbdelem[bid] ; id++ ) {


	    // Node id

	    uint node = mesh.bd.bdPoints[bid][id];

		
	    // Move over velocities and check for non-existing neighbour

	    for( uint k = 0 ; k < q ; k++ ) {

		// This velocity sets virtual node position
		// points[node] + vel[k]
		
		    
		if( mesh.nb[node][ reverse[k] ] == -1 ) {


		    // First compute distance to other nodes conecting "node"
		    // Virtual   node: node + v
		    // Neighbour node: node + v'
		    // Possible distance between virtual and neighbours = (node + v) - (node + v') = v - v'
		    // In order to avoid problems with periodic boundaries, mesh points are used to compute distance
		    //
		    // Velocity sets may not fill all closest neighbours (as in D3Q15), so also check for neighbours of neighbours


		    // Distance array

		    vector< vector<int> > distance;
		    

		    // Find closest nodes from vtkCells

		    for( auto cell : mesh.nodeToCells[node] ) {

			for( auto nv : mesh.vtkCells[cell] ) {

			    // Separation vector

			    int sep[3] = {  (mesh.points[nv][0]  - mesh.points[node][0] - vel[k][0]),
				  	    (mesh.points[nv][1]  - mesh.points[node][1] - vel[k][1]),
					    (mesh.points[nv][2]  - mesh.points[node][2] - vel[k][2])  };
					 

			    // Squared distance
			    
			    int dist = sep[0]*sep[0] + sep[1]*sep[1] + sep[2]*sep[2];

			    if( dist <= 3 ) {


				// Find node along sep using cell info
				
				int inNode = nv;

				for( auto nvCell : mesh.nodeToCells[nv] ) {

				    for( auto vrt : mesh.vtkCells[nvCell] ) {

					int otherSep[3] = {  (mesh.points[vrt][0]  - mesh.points[nv][0]),
							     (mesh.points[vrt][1]  - mesh.points[nv][1]),
							     (mesh.points[vrt][2]  - mesh.points[nv][2])  };

					if( otherSep[0] == sep[0] ) {
					    
					    if( otherSep[1] == sep[1] ) {

						if( otherSep[2] == sep[2] ) {

						    inNode = vrt;

						}

					    }

					}

				    }

				}
				


				// Append to distance array
				
				distance.push_back( {dist, nv, inNode} );

			    }
			    
			}

		    }
		    



		    // Sort distance

		    // std::sort( distance.begin(), distance.end(), compareDist );
		    for( uint i = 0 ; i < distance.size() ; i++ ) {

			if( distance[i][0] < distance[0][0] ) {

			    for( uint j = 0 ; j < 3 ; j++ )
				distance[0][j] = distance[i][j];

			}

		    }
		    


		    // Add to nodes

		    vnodes.push_back( { (int)node, (int)k, distance[0][1], distance[0][2]} );
						

		}

	    }
		

	}
	 
	   	    

    }
  

    
    
}
