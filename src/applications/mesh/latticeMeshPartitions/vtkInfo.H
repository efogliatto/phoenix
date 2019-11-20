#ifndef VTKINFO_H
#define VTKINFO_H

/**
 * @file vtkInfo.h
 * @author Ezequiel O. Fogliatto
 * @date 15 Dic 2017
 * @brief VTK output information
 */

/**
 * @brief VTK output information
 *
 * Information related to VTK uotput
 */

#include <dataTypes.h>

#include <vector>

#include <string>

typedef struct {

    /** Total number of scalar fields */
    uint nscalar;   

    /** Total number of vector fields */
    uint nvector;   

    /** Total number of pdf fields */
    uint npdf;   

    /** Scalar fields names*/
    std::vector< std::vector<std::string> > scalarFields;   

    /** Vector fields names*/
    std::vector< std::vector<std::string> > vectorFields;

    /** Pdf fields names*/
    std::vector< std::vector<std::string> > pdfFields;
    
} vtkInfo;

#endif // VTKINFO_H
