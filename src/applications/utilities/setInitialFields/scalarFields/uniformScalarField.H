#ifndef UNIFORMSCALARFIELD_H
#define UNIFORMSCALARFIELD_H

#include <scalarField.H>

#include <latticeMesh.H>

void uniformScalarField( scalarField& field, const latticeMesh& mesh, const std::string& fname, const std::string& sptype );

#endif // UNIFORMSCALARFIELD_H
