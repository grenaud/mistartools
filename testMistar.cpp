/*
 * testMistar.cpp
 * Date: Jan-25-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "MistarParser.h"

using namespace std;

int main (int argc, char *argv[]) {



    string line;
    igzstream myFile;
    string filename = string(argv[1]);
    myFile.open(filename.c_str(), ios::in);
    vector<string> lines;
    vector<string> populationNames;
    unsigned int numberPopulations;
   

   
    if (myFile.good()){

	while ( getline (myFile,line)){
	    //	    cout<<line<<endl;
	    if(line[0] == '#'){
	       
		if(strBeginsWith(line, "#chr")){
		    vector<string> fields=allTokens(line,'\t');
		    if(fields[0] != "#chr")   { cerr<<"Field #1 of header must be #chr ";    exit(1); }
		    if(fields[1] != "coord")  { cerr<<"Field #2 of header must be coord ";   exit(1); }
		    if(fields[2] != "REF,ALT"){ cerr<<"Field #3 of header must be REF,ALT "; exit(1); }
		    if(fields[3] != "root")   { cerr<<"Field #4 of header must be root ";    exit(1); }
		    if(fields[4] != "anc")    { cerr<<"Field #5 of header must be anc ";     exit(1); }

		    for(unsigned int i=3;i<fields.size();i++){
			populationNames.push_back(fields[i]);
			numberPopulations++;
		    }

		    //break;
		}else{
		  
		}
	    
	    }else{

		lines.push_back(line);
	    }
	}

	myFile.close();
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	return 1;
    }
    cout<<"done first read\t"<<lines.size()<<endl;
    AlleleRecords * test;

    MistarParser mp1 (&lines,populationNames);
    cout<<"done constr"<<endl;
    while(mp1.hasData()){
	cout<<"data"<<endl;
	test = mp1.getData();
	cout<<test->coordinate<<endl;
    }
    cout<<"done read"<<endl;
    return 0;



    
    cout<<"opening "<<argv[1]<<endl;
    
    MistarParser mp (argv[1]);

    cout<<"Header\n\n\n\n "<<mp.getHeader()<<endl;
    cout<<"\n\n\n\n "<<endl;
    while(mp.hasData()){
	cout<<"data"<<endl;
	test = mp.getData();
	cout<<test->coordinate<<endl;
    }


    return 0;
}

