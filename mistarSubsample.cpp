/*
 * mistaruniq
 * Date: Apr-23-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

// #define DEBUG

#include "utils.h"
#include "MistarParser.h"
//#include "mistarunion.h"
#include "mistarOperations.h"

using namespace std;



int main (int argc, char *argv[]) {
    

    if(argc < 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr<<"usage: "<<argv[0]<<" [fraction] [mistar file 1] ...\nwill print the line of the original mistar file with probability [fraction]"<<endl;
	return 1;       
    }

    



    
    MistarParser  mp        (argv[argc-1]);
    string fraction = string(argv[argc-2]);
    double fractionD = destringify<double>(fraction);
    if(fractionD<0 || fractionD>1){
	cerr<<"fraction: "<<argv[argc-2]<<" must be between 0 and 1"<<endl;
	return 1;	
    }

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    cout<<"#MISTAR"<<endl;    	
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;
    cout<<"#MISTARSUBSAMPLE:"<<endl;

    cout<<"#MISTARSUBSAMPLE#1"<<endl;
    cout<<""<<mp.getHeader("#\t")<<"\n";
    cout<<"#chr\tcoord\tREF,ALT"<<"\t"<<vectorToString(*mp.getPopulationsNames(),"\t")<<endl;


    AlleleRecords * dataRow;

    while(mp.hasData()){
	dataRow = mp.getData();
	if(randomProb()<=fractionD)
	    cout<<*dataRow<<endl;	    
    }




    cerr<<"Program finished gracefully"<<endl;



    return 0;
}

