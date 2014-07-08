/*
 * mistarfilter
 * Date: Mar-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>

#include "utils.h"

#include "MistarParser.h"
#include "GenomicRange.h"

using namespace std;


int main (int argc, char *argv[]) {

    string usage=string(""+string(argv[0])+"  [mistar file]"+
			"\nThis program will print the records in the mistar file as bed file\n\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }
    
    

    MistarParser mp   (argv[argc-1]);
    string chrName="-1";
    unsigned int previousCoordinate=0;
    unsigned int startBed=0;
    unsigned int endBed=0;
    

    AlleleRecords * dataRow;

	    
    while(mp.hasData()){
	dataRow = mp.getData();
	// cout<<"test\t"<<(*dataRow)<<endl;
	
	if(dataRow->chr != chrName){//new chr

	    if(chrName != "-1"){
		cout<<chrName<<"\t"<<startBed<<"\t"<<endBed<<endl;
	    }

	    previousCoordinate = dataRow->coordinate;
	    chrName            = dataRow->chr;
	    startBed = dataRow->coordinate-1;
	    endBed   = dataRow->coordinate;


	}else{
	 
	    if(previousCoordinate >= dataRow->coordinate){
		cerr<<"There seems to be a unsorted coordinate in the mistar file, needs to be sorted coordinate: "<<previousCoordinate<<endl;
		return 1;
	    }

	    if( (previousCoordinate+1) == dataRow->coordinate){
		endBed++;
	    }else{
		cout<<chrName<<"\t"<<startBed<<"\t"<<endBed<<endl;
		startBed= dataRow->coordinate-1;
		endBed  = dataRow->coordinate;
	    }

	    previousCoordinate = dataRow->coordinate;

	}
    }	
    cout<<chrName<<"\t"<<startBed<<"\t"<<endBed<<endl;

	
    return 0;
}

