#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <string>

#include <iostream>

#include <fstream>

#include <vector>

#include <dataTypes.H>



class dictionary {

protected:

    
    /** File content */

    std::vector< std::string > content;

    
    /** Braced entry */
    
    const std::vector<std::string> bracedEntry( const std::string& ename, const std::vector<std::string>& elist ) const;

    
    /** Tokenized entry */

    const std::vector<std::string> tkEntry( const std::string& ename ) const;

    
    /** Extract string from entry */

    const void extract( const std::vector<std::string>& elist, std::string& res ) const;


    /** Extract scalar from entry */

    const void extract( const std::vector<std::string>& elist, scalar& res ) const;


    /** Extract int from entry */

    const void extract( const std::vector<std::string>& elist, int& res ) const;


    /** Extract scalar vector from entry */

    const void extract( const std::vector<std::string>& elist, std::vector<scalar>& res ) const;


    /** Extract vector int from entry */

    const void extract( const std::vector<std::string>& elist, std::vector<int>& res ) const;    

    


public:


    // Constructors and destructors

    /** Default constructor */

    dictionary(const std::string& fname);


    /** Default destructor */

    ~dictionary();




    // Look single entry
    template<typename T>
    const T lookUp(const std::string& entry) const {

	T out;

	std::vector<std::string> elist = tkEntry(entry);

	if( elist.size() == 0 ) {

	    std::cout << " [ERROR]  Unable to find entry " << entry << std::endl;

	    exit(1);

	}

	else {

	    extract( elist, out );

	}
	
	return out;

    }

};

#endif // DICTIONARY_H
