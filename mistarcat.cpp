/*
 * mistarcat
 * Date: Feb-04-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>


#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {

    if(argc == 1){
	cerr << "This program concatenates many files where the header is found in the\nfirst file and does not use the headers from the remaining ones\nIt prints to the /dev/stdout\nUsage: "<<argv[0]<<" [mistar file#1] [mistar file#2] ... "<<endl;
       return 1;	
    }



    for(int i=1;i<argc;i++){

	string line;
	igzstream myFile;
	string filename = string(argv[i]);
	//cerr<<filename<<endl;
	myFile.open(filename.c_str(), ios::in);
	bool header=true;
	bool haveseenHeader=false;

	if (myFile.good()){
	    while ( getline (myFile,line)){
		//cerr<<line<<endl;
		//		if(i!=1){
		
		if(line[0] == '#'){
		    header=true;

		    if(haveseenHeader){
			cerr << "ERROR: Line "<<line<<" in "<<filename<<" cannot start with a #"<<endl;
			return 1;
		    }
		}
		
		if(header){
		    if(i==1){ //use header from first file
			cout<<line<<endl;
		    } //else do nothing
		    
		    header=false;
		}else{		    
		    haveseenHeader=true;
		    cout<<line<<endl;
		}
	    }
	    myFile.close();
	}else{
	    cerr << "Unable to open file "<<filename<<endl;
	    return 1;
	}
    }
    return 0;
}

