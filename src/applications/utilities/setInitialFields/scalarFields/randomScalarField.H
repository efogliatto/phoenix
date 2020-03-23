#ifndef RANDOMSCALARFIELD_H
#define RANDOMSCALARFIELD_H

#include <scalarField.H>

#include <latticeMesh.H>

void randomScalarField( scalarField& field, const latticeMesh& mesh, const std::string& fname, const std::string& sptype );

#endif // RANDOMSCALARFIELD_H
