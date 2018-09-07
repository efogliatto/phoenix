#ifndef TIMEOPTIONS_H
#define TIMEOPTIONS_H

#include <dataTypes.H>

#include <string>

#include <dictionary.H>

#include <chrono>

#include <tuple>



class timeOptions {


protected:


    /** Processor id */

    uint pid;
    
    
    /** Start time */
    
    uint start;   

    
    /** End time */
    
    uint end;   


    /** Current time */
    
    uint current;   


    /** Write interval (write every writeInterval steps) */
    
    uint writeInterval;   


    
    /** Start time for elapsed time */
    
    typedef std::chrono::high_resolution_clock clock_;

    typedef std::chrono::duration<double, std::ratio<1> > second_;

    std::chrono::time_point<clock_> beg;
    


    /** Time step inside write interval */
    
    uint stp;   


    /** Data format */

    std::string dataFormat;



    /** Update ensight case file */

    const void updateCaseFile() const;


    /** Simulation files */

    std::vector<std::string> scalarFields;

    std::vector<std::string> vectorFields;

    std::vector< std::tuple<std::string, uint> > pdfFields;    




public:

    
    /** Default constructor */

    timeOptions( uint id = 0 );


    /** Default destructor */

    ~timeOptions();



    /** Seconds since class instantiation*/

    const scalar elapsed() const;


    /** Update time: move to next time step */

    const bool update();


    /** Write flag */

    const bool write() const;


    /** Add simulation file names */

    const void addScalarField( const std::string& name );

    const void addVectorField( const std::string& name );

    const void addPdfField( const std::string& name, const uint& q = 9 );    

};

#endif // TIMEOPTIONS_H
