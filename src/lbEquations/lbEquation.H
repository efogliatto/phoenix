#ifndef LBEQUATION_H
#define LBEQUATION_H

#include <scalarField.H>

#include <vectorField.H>

#include <pdfField.H>

#include <map>


class lbEquation {

protected:

    
    /** Equation model name */

    std::string ename;


    /** Reference to lattice mesh */

    const latticeMesh& mesh;


    /** Reference to time options */

    timeOptions& Time;


    /** Cached lattice q (efficiency) */

    const uint _q;
    

    /** Reference to external pdf field */

    pdfField& pdf;


    /** Swap auxiliary field */

    std::vector< std::vector<scalar> > swap;    
    



public:


    /** Default constructor */

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