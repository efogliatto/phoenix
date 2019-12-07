#ifndef UPDATEPOINTSANDCELLS_H
#define UPDATEPOINTSANDCELLS_H

#include <dataTypes.H>

#include <vector>

void updatePointsAndCells( std::vector< std::vector<uint> >& basePoints,
			   std::vector< std::vector<uint> >& baseCells,
			   std::vector< std::vector<uint> >& newPoints,
			   std::vector< std::vector<uint> >& newCells,
			   const std::vector<bool>& isInside );

#endif // UPDATEPOINTSANDCELLS_H
