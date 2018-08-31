#ifndef DATATYPES_H
#define DATATYPES_H


/**
 * @file dataTypes.h
 * @author Ezequiel O. Fogliatto
 * @date 31 Aug 2018
 * Basic pre-defined constants and data types
 */

/** Unsigned int short name */

typedef unsigned int uint;


/** Scalar for precision change */

#ifdef DP

typedef double scalar;

#elif SP

typedef float scalar;

#endif



#endif // DATATYPES_H
