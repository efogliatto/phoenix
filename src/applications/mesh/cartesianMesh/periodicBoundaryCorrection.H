#ifndef PERIODICBOUNDARYCORRECTION_H
#define PERIODICBOUNDARYCORRECTION_H

#include <dataTypes.H>

#include <vector>

#include <unordered_map>

#include <string>

#include <utility>

#include <periodicBnds.H>


void periodicBoundaryCorrection( std::vector< std::vector<int> >& nb,
				 const std::vector<periodicBnds>& periodicPairs,
				 const std::vector< std::vector<uint> >& points,
				 const std::unordered_map< std::string, std::vector<uint> >& boundaries );

#endif // PERIODICBOUNDARYCORRECTION_H
