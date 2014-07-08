/*
 * mistarsplit
 * Date: Apr-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "MistarParser.h"

using namespace std;

int main (int argc, char *argv[]) {
    string usage=string(""+string(argv[0])+" <options> [mistar file] [outprefix] "+
			"\n\nThis program takes a mistar matrix and split it according to chromosome\n\n"+
			"\tOptions\n"+
			// "// \t\t"+"--chrcoord"+"\t"+"Print the chr/coord (Default: "+stringify(printChr)
			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }



    //starts at 1 and except the last two
    for(int i=1;i<(argc-2);i++){ 
	// if(string(argv[i]) == "--chrcoord" ) {
	//     printChr=true;
	//     continue;
	// }

        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1;           

    }

    unsigned int totalRecords=0;
    unsigned int keptRecords=0;

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    // cout<<"#MISTAR"<<endl;
    // cout<<"#PG:"<<programLine<<endl;
    // cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    // cout<<"#DATE: "<<getDateString()<<endl;

    MistarParser mp   (argv[argc-2]);
    string commonPrefix = string(argv[argc-2]);

    AlleleRecords * dataRow;
    string lastChr="##";

    while(mp.hasData()){
	dataRow = mp.getData();
	if(dataRow->chr != lastChr){
	    
	    lastChr = dataRow->chr ;
	}else{
	    
	}
    }


    return 0;
}

