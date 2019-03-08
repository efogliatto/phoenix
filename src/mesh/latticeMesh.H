#ifndef LATTICEMESH_H
#define LATTICEMESH_H

#include <vector>

#include <mpiInfo.H>

#include <latticeModelCreator.H>

#include <map>


class latticeMesh {

    
protected:


    /** Total number of points per patch */

    uint nPoints;


    /** Points coordinates */

    std::vector< std::vector<int> > points;


    /** Neighbour indices array */

    std::vector< std::vector<int> > nb;
    

    /** MPI comm information */

    const mpiInfo parallel;

    
    /** Lattice Model */

    latticeModel* lbmodel;


    /** Boundary nodes */

    std::map< std::string, std::vector<uint> > boundary;


    /** Tell if node is on boundary. Performance */

    std::vector<uint> isOnBnd;


    /** Boundary to which node belongs */

    std::map<uint, std::string> nodeToBnd;


    /** Aplicable nodes */
    
    std::map< std::string, std::vector<uint> > closestNodes;
    



    

    


    /** Read mesh points */

    const void readPoints();


    /** Read neighbours */

    const void readNeighbours();

    
    /** Read boundary nodes */

    const void readBoundaryNodes();
    



public:

    
    /** Default constructor */

    latticeMesh( const int& pid, const bool& msg = true );


    /** Default destructor */

    ~latticeMesh();


    /** Reference to mpiInfo */

    const mpiInfo& mpi() const { return parallel; }

    
    /** Process id */

    const int pid() const { return parallel.id(); }

    
    /** Total number of processors */

    const int wsize() const { return parallel.wsize(); }
    

    /** Total number of points */

    const uint npoints() const { return nPoints; }


    /** Number of local points */

    const uint local() const { return parallel.local(); }       
    
    
    /** Nodes per patch */

    const uint npp( const uint& id ) const { return parallel.npp(id); }

    
    /** Reference to model */

    const latticeModel* lmodel() const { return lbmodel; }


    /** Reference to neighbours array */

    const std::vector< std::vector<int> >& nbArray() const { return nb; }


    /** Reference to boundary indices */

    const std::vector<uint>& boundaryNodes( const std::string& bdname ) const;

    const std::map< std::string, std::vector<uint> >& boundaries() const { return boundary; }


    /** Lattice point coordinates */

    const std::vector<int>& latticePoint( const uint& id ) const { return points[id]; }


    /** Check if point is on boundary */

    const uint isOnBoundary( const uint id ) const { return isOnBnd[id]; };


    /** Boundary which node belongs */

    const std::string nodeToBoundary( const uint id ) const;


    /** Nodes close to boundaries */

    const std::map< std::string, std::vector<uint> >& nodesCloseToBoundary() const { return closestNodes; }


    /** Boundary normal type */

    enum class normalType{ X0, X1, Y0, Y1, Z0, Z1, UNDEF };    
    
    
};

#endif // LATTICEMESH_H
