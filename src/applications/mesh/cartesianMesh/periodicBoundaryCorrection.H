#ifndef PERIODICBOUNDARYCORRECTION_H
#define PERIODICBOUNDARYCORRECTION_H

#include <dataTypes.H>

#include <vector>

#include <map>

#include <string>

#include <utility>


void periodicBoundaryCorrection( std::vector< std::vector<int> >& nb,
				 const std::map< std::pair<std::string,std::string>, std::vector<scalar> >& periodicPairs,
				 const std::vector< std::vector<uint> >& points,
				 const std::map< std::string, std::vector<uint> >& boundaries );

#endif // PERIODICBOUNDARYCORRECTION_H
