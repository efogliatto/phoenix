#include <findClosestBoundary.H>

#include <CGAL/AABB_tree.h>

#include <CGAL/AABB_traits.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/algorithm.h>




typedef Kernel::Point_3 Point;

typedef Kernel::FT FT;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;

typedef CGAL::AABB_traits<Kernel, Primitive> Traits;

typedef CGAL::AABB_tree<Traits> Tree;





using namespace std;

void findClosestBoundary( unordered_map< string, vector<uint> >& boundaries,
			  const vector< vector<uint> >& meshPoints,
			  const vector< vector<int> >& nb,
			  const vector< pair<string, Polyhedron> >& bdPolyMap ) {



    // Distance map

    map< uint, unordered_map<string, FT> > dstMap;

    
    // Priority weight

    map<string, uint> priority;

    for( uint i = 0 ; i < bdPolyMap.size() ; i++ )
	priority[ bdPolyMap[i].first ] = i;
   

    
    // Constants

    const uint nPoints = meshPoints.size();

    const uint q = nb[0].size();
    

    // Move over neighbours

    for( uint i = 0 ; i < nPoints ; i++ ) {

    	bool isOnBnd(false);

    	for( uint k = 0 ; k < q ; k++ ) {

    	    if ( nb[i][k] == -1 )
    		isOnBnd = true;

    	}


	if( isOnBnd ) {

	    for( auto poly : bdPolyMap )
		dstMap[i][poly.first] = 0;


	}

    }


    

    // Construct AABB trees for each boundary

    for( auto& poly : bdPolyMap ) {

    	Tree tree( faces(poly.second).first, faces(poly.second).second, poly.second );

    	tree.accelerate_distance_queries();


    	for( auto& dst : dstMap ) {

    	    Point query( meshPoints[dst.first][0], meshPoints[dst.first][1], meshPoints[dst.first][2] );

    	    FT sqd = tree.squared_distance(query);

    	    dst.second[ poly.first ] = tree.squared_distance(query);

    	}	

    }
  


    // Check closest surface

    for( auto& dst : dstMap ) {

	FT minDist = 1000000;

	string closest;

	for( auto& surf : dst.second ) {

	    if( surf.second < minDist ) {

		minDist = surf.second;

		closest = surf.first;

	    }

	    else {

		if( surf.second == minDist ) {

		    if( priority[surf.first] > priority[closest] ) {		    

			minDist = surf.second;

			closest = surf.first;

		    }

		}

	    }

	}


	// Need to sort according to priority

	


	// Append to boundaries dictionary

	boundaries[closest].push_back( dst.first );

    }

    
    

}
