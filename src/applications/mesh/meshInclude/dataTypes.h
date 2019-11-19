#ifndef DATATYPES_H
#define DATATYPES_H


/**
 * @file dataTypes.h
 * @author Ezequiel O. Fogliatto
 * @date 26 Apr 2018
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


/** Field creation options */

typedef enum {

    MUST_READ,
    
    NO_READ    

} fieldOpt;


/** Data format */

typedef enum {

    asciiRaw,

    binaryRaw,

    pvtu,

    ensight

} dataFormat;


#endif // DATATYPES_H
