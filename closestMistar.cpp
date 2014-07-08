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

    if(argc == 1 ||
       (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
       ){
	cerr<<"usage: "<<argv[0]<<" <options> [mistar file]\nwill print to stdout the distance to the closest site for each record"
	    <<endl;
	return 1;
    }
    


    MistarParser mp (argv[argc-1]);
    AlleleRecords * test;
    unsigned int totalRecords=0;
    unsigned int lastCoordinateM2=0;
    unsigned int lastCoordinateM1=0;
    unsigned int lastCoordinateM =0;
    string chr;

    // if(mp.hasData()){
    // 	test = mp.getData();
    // 	lastCoordinateM2=test->coordinate;
    // 	chr=test->chr;
    // 	totalRecords++;
    // }else{
    // 	cerr<<"Program "<<argv[0]<<" looked at "<<totalRecords<<" terminated gracefully"<<endl;
    // 	return 0;
    // }

    // if(mp.hasData()){
    // 	test = mp.getData();
    // 	lastCoordinateM1=test->coordinate;
    // 	if( chr != test->chr){
    // 	    cerr<<"Chromosomes differ for line :"<<(*test)<<endl;
    // 	    return 1;
    // 	}
    // 	totalRecords++;
    // }else{
    // 	cerr<<"Program "<<argv[0]<<" looked at "<<totalRecords<<" terminated gracefully"<<endl;
    // 	return 0;
    // }

    // cout<<(lastCoordinateM1-lastCoordinateM2)<<endl;

    bool newChr1=true;
    bool newChr2=true;


    while(mp.hasData()){
	//cout<<"data"<<endl;
	test = mp.getData();

	if( chr != test->chr){
	    //cerr<<"Chromosomes differ for line :"<<(*test)<<endl;
	    newChr1=true;
	    newChr2=true;
	    lastCoordinateM2=0;
	    lastCoordinateM1=0;
	    lastCoordinateM =0;
	    //return 1;
	}

	if(newChr1){
	    lastCoordinateM2=test->coordinate;
	    chr=test->chr;
	    totalRecords++;
	    newChr1=false;

	    continue;
	}

	if(newChr2){
	    lastCoordinateM1=test->coordinate;
	    if(lastCoordinateM1<lastCoordinateM2){
		cerr<<"Coordinate are not sorted :"<<(*test)<<endl;
		return 1;
	    }

	    chr=test->chr;
	    totalRecords++;
	    cout<<(lastCoordinateM1-lastCoordinateM2)<<endl;
	    newChr2=false;

	    continue;
	}


	lastCoordinateM=test->coordinate;

	if(lastCoordinateM<lastCoordinateM1){
	    cerr<<"Coordinate are not sorted :"<<(*test)<<endl;
	    return 1;
	}

	if( (lastCoordinateM1-lastCoordinateM2) < (lastCoordinateM-lastCoordinateM1) ){
	    cout<<(lastCoordinateM1-lastCoordinateM2)<<endl;
	}else{
	    cout<<(lastCoordinateM-lastCoordinateM1)<<endl;
	}
	
	//for next iteration
	lastCoordinateM2=lastCoordinateM1;
	lastCoordinateM1=lastCoordinateM;

	totalRecords++;
    }
    cout<<(lastCoordinateM1-lastCoordinateM2)<<endl;
    
    cerr<<"Program "<<argv[0]<<" looked at "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

