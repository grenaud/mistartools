/*
 * mistar2binary
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
    string usage=string(""+string(argv[0])+" <options> [mistar file] "+
			"\n\nThis program takes a mistar matrix and prints some stats\n\n"+
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

    MistarParser mp   (argv[argc-1]);
    AlleleRecords * dataRow;

    while(mp.hasData()){
	dataRow = mp.getData();
	totalRecords++;
    }

    cerr<<"Program "<<argv[0]<<" saw "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully";

    return 0;
}

