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

    const bool write(const bool updCaseFile=true) const;


    /** Add simulation file names */

    const void addScalarField( const std::string& name );

    const void addVectorField( const std::string& name );

    const void addPdfField( const std::string& name, const uint& q = 9 );


    /** Start time */

    const uint startTime() const { return start; }


    /** Current time */

    const uint currentTime() const { return current; }    


    /** Match time to ensight index */

    const uint timeToIndex( const uint& tid ) const;

    
    /** Update ensight case file */

    const void updateCaseFile() const;


    /** Time-index tuple vector */

    const std::vector<uint> timeList() const;


    /** Check lbm.case and save other fields that are not in register */

    void keepRegisteredFields();    
    
};




// Read - write options

enum class IO : int {MUST_READ, NO_READ, MUST_WRITE, NO_WRITE};



#endif // TIMEOPTIONS_H
