#ifndef COMPUTENEIGHBOURSFROMCELLS_H
#define COMPUTENEIGHBOURSFROMCELLS_H

#include <dataTypes.H>

#include <vector>

#include <latticeModelCreator.H>


void computeNeighboursFromCells( std::vector< std::vector<int> >& nb, std::vector< std::vector<uint> >& points, std::vector< std::vector<uint> >& cells, const latticeModel* lbmodel );

#endif // COMPUTENEIGHBOURSFROMCELLS_H
