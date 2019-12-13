#ifndef SYMMETRYPLANESCORRECTION_H
#define SYMMETRYPLANESCORRECTION_H

#include <dataTypes.H>

#include <vector>

#include <unordered_map>

#include <string>

#include <latticeModelCreator.H>

void symmetryPlanesCorrection( std::vector< std::vector<int> >& nb,
			       const std::string& plane,
			       const scalar offset,			       
			       const std::vector< std::vector<uint> >& points,
			       const std::unordered_map< std::string, std::vector<uint> >& boundaries,
			       const latticeModel* lbmodel );

#endif // SYMMETRYPLANESCORRECTION_H
