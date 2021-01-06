#ifndef TGPROFILESCALARFIELD_H
#define TGPROFILESCALARFIELD_H

#include <scalarField.H>

#include <latticeMesh.H>

void tgProfileScalarField( scalarField& field, const latticeMesh& mesh, const std::string& fname, const std::string& sptype );

#endif // TGPROFILESCALARFIELD_H
