#include <dictionary.H>
#include <sstream>
#include <boost/algorithm/string/replace.hpp>



using namespace std;



/** Default constructor */

dictionary::dictionary(const std::string& fname) {


    // Open file and load content

    ifstream inFile;

    inFile.open( fname.c_str() );

    if( inFile.is_open() == false ){
	
    	cout << " [ERROR]  Unable to find file " << fname << endl;
	
    	exit(1);
	
    }
    

    // Read line by line
    
    string line;

    vector<string> tmp;
    

    while( !inFile.eof() ) {

    	// Read line
	
    	getline(inFile, line);


	// Replace specific tokens

	boost::replace_all(line, "(", " ( ");
	boost::replace_all(line, ")", " ) ");
	boost::replace_all(line, "{", " { ");
	boost::replace_all(line, "}", " } ");
	boost::replace_all(line, "/*", " /* ");
	boost::replace_all(line, "*/", " */ ");
	boost::replace_all(line, ";", " ; ");
	
	

	// Split line into words and tokens, and remove whitespaces
	    
	istringstream ss(line);

	do {
		
	    string sub;

	    ss >> sub;
		
	    if(!sub.empty())	
		tmp.push_back( sub );
		    
		
	} while(ss);
	    
	
    }


    // Close file
    
    inFile.close();



    // Remove comments

    for( uint i = 0 ; i < tmp.size()-1 ; i++ ) {

	if( tmp[i] == "/*" ) {

	    do {

		i++;

	    } while (tmp[i] != "*/");

	    i++;
	    
	}

	content.push_back( tmp[i] );

    }
    



}


/** Default destructor */

dictionary::~dictionary() {}





/** Braced entry */
    
const vector<string> dictionary::bracedEntry( const string& ename, const vector<string>& elist ) const {

    vector<string> out;


    for( uint i = 0 ; i < elist.size()-1 ; i++ ) {

	if(  ( elist[i] == ename )   &&  ( elist[i+1] == "{" ) ) {

	    i = i+2;
	    
	    uint count = 1;

	    while(  (count != 0)  &&  ( i < elist.size() )  ) {

	    	if(elist[i] == "{")
		    count++; 

		if(elist[i] == "}")
		    count--;

		if(count != 0)
		    out.push_back( elist[i] );

		i++;
		
	    }	   

	}

    }
    

    return out;

}



/** Tokenized entry */

const vector<string> dictionary::tkEntry( const std::string& ename ) const {


    // Tokenize entry

    vector<string> tokens;
    
    string token;

    istringstream tokenStream(ename);
    
    while ( getline(tokenStream, token, '/') )    {
	
      tokens.push_back(token);
      
    }



    // Look until last entry

    vector<string> out( content );
    
    for(uint i = 0 ; i < tokens.size()-1 ; i++ ) {

	out = bracedEntry( tokens[i], out );

    }



    vector<string> res;

    for(uint i = 0 ; i < out.size() ; i++) {

	if( out[i] == tokens.back() ) {

	    i++;
	    
	    while( out[i] != ";" ) {

		res.push_back( out[i] );

		i++;

	    }

	}

    }


    return res;

}




/** Extract string from entry */

const void dictionary::extract( const std::vector<std::string>& elist, string& res ) const {

    res = elist[0];

}



/** Extract scalar from entry */

const void dictionary::extract( const std::vector<std::string>& elist, scalar& res ) const {

    res = (scalar)stod( elist[0] );

}


/** Extract int from entry */

const void dictionary::extract( const std::vector<std::string>& elist, int& res ) const {

    res = stoi( elist[0] );

}



/** Extract scalar vector from entry */

const void dictionary::extract( const std::vector<std::string>& elist, vector<scalar>& res ) const {

    vector<scalar> out;
    
    for( uint i = 0 ; i < elist.size() ; i++ ) {

    	if( elist[i] == "(" ) {

    	    i++;
	    
    	    while( elist[i] != ")" ) {

    		out.push_back( stod(elist[i]) );

		i++;

    	    }

    	    i = elist.size();

    	}

    }


    res = out;


}
