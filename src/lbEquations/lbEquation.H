#ifndef LBEQUATION_H
#define LBEQUATION_H

#include <scalarField.H>

#include <vectorField.H>

#include <pdfField.H>

#include <map>

#include <relaxModelCreator.H>


/**
 * @file lbEquation.H
 * @author Ezequiel O. Fogliatto
 * @date 07 Jan 2020
 * @brief Base class for LBE
 */

/**
 *  A test class. A more elaborate class description.
 */

class lbEquation {

protected:
    

    std::string ename;           /**< Equation model name. */

    const latticeMesh& mesh;     /**< Reference to lattice mesh */
    
    timeOptions& Time;           /**< Reference to time options */
     
    pdfField& _pdf;              /**< Reference to external pdf field */  

    std::vector<scalar> _Tau;    /**< Relaxation coefficients */
        
    std::vector< std::vector<scalar> > _swap;  /**< Swap auxiliary field */

    std::unique_ptr<relaxModel> _relax;      /**< Relaxation model */
    



public:


    /** 
     * Default constructor. 
     * @param name equation name.
     * @param mesh_ external reference to mesh object.
     * @param Time_ external reference to time object.
     * @param pdf_ external reference to distribution field.
     */
    lbEquation( const std::string& name,
		const latticeMesh& mesh_,
		timeOptions& Time_,
		pdfField& pdf_);


    /** Default destructor */

    ~lbEquation();



    /** Collision process */

    virtual const void collision() = 0;


    /** Streamming process */

    const void streaming();


    /** Set pdf to equilibrium values */

    virtual const void setEquilibrium() = 0;
    

};

#endif // LBEQUATION_H
