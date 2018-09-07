#ifndef TIMEOPTIONS_H
#define TIMEOPTIONS_H

#include <dataTypes.H>

#include <string>

#include <dictionary.H>

#include <chrono>



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

};

#endif // TIMEOPTIONS_H
