/*
 * mistar2segtor
 * Date: Aug-16-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "MistarParser.h"

using namespace std;

int main (int argc, char *argv[]) {
    // bool printChr=false;

    string usage=string(""+string(argv[0])+"  [mistar file] "+
			"\n\nThis program takes a mistar matrix and prints input for Segtor\n\n");
			// "\tOptions\n"+
			// "\t\t"+"--chrcoord"+"\t"+"Print the chr/coord (Default: "+stringify(printChr)+" )\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }



    //starts at 1 and except the last two
    for(int i=1;i<(argc-1);i++){ 
	// if(string(argv[i]) == "--chrcoord" ) {
	//     printChr=true;
	//     continue;
	// }

        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1;           

    }

    unsigned int totalRecords=0;
    // unsigned int keptRecords=0;

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }


    MistarParser mp   (argv[argc-1]);

    AlleleRecords * dataRow;

    while(mp.hasData()){
	dataRow = mp.getData();
	
	cout<<dataRow->chr<<"\t"<<dataRow->coordinate<<"\t"<<dataRow->ref<<"\t"<<dataRow->alt<<"\tsnp"<<totalRecords<<endl;
	totalRecords++;
    }

    cerr<<"Program "<<argv[0]<<" wrote "<<totalRecords<<" terminated gracefully";

    return 0;
}

