/*
 * usepopasrootanc.cpp
 * Date: Feb-22-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <set>

#include "utils.h"
#include "MistarParser.h"

using namespace std;

int main (int argc, char *argv[]) {

    if(argc == 1 ||
       (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
       ){
	cerr<<"usage: "<<argv[0]<<" <options> [mistar file] \"poproot\" \"popanc\"\nwill print to stdout"
	    <<endl
	    <<"This program will use specified populations as root and ancestor and produce lines with only those two populations"<<endl
	    <<"ex: "<<argv[0]<<" data.mst rootpop ancpop"<<endl
	    <<"Options:"<<endl
	    // <<"\t--justtransv\t\tOnly allow transversions (Default "+boolStringify(limitToTransversions)+" )\n"
	    // <<"\t--noprivate\t\tDo not allow private mutations (Default "+boolStringify(noprivate)+" )\n"

	    <<endl;
	
	return 1;
    }
    
    //all but two last
    for(int i=1;i<(argc-3);i++){ 
	

        // if( string(argv[i]) == "--justtransv"  ){
	//     limitToTransversions=true;
        //     continue;
        // }

        // if( string(argv[i]) == "--noprivate"  ){
	//     noprivate=true;
        //     continue;
        // }

	cerr<<"Error unknown option "<<argv[i]<<endl;
	exit(1);
    }

    MistarParser mp                  (string(argv[argc-3]));
    string poproot           =string(argv[argc-2]);
    string popanc            =string(argv[argc-1]);
  

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    cout<<"#MISTAR"<<endl;    	
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;
    cout<<"#USEPOPASROOTANC:"<<endl;
    //for(unsigned int i=0;i<vectorOfMP.size();i++){ 
    cout<<"#USEPOPASROOTANC#"<<endl;
    cout<<""<<mp.getHeader("#\t")<<"\n";
    //}
    
    cout<<"#chr\tcoord\tREF,ALT\troot\tanc"<<endl;

    //    vector<string> newPopVector;
    unsigned int rootidx =0;
    unsigned int ancidx  =1;
    bool rootflg = false;
    bool ancflg  = false;

    for(unsigned int i=0;i<mp.getPopulationsNames()->size();i++){

	if(mp.getPopulationsNames()->at(i) ==  poproot){
	    rootidx =i;
	    rootflg=true;
	}

	if(mp.getPopulationsNames()->at(i) ==  popanc){
	    ancidx =i;
	    ancflg=true;
	}

    }

    if(!rootflg ){
	cerr<<"Error: the root pop was not found"<<endl;
	return 1;
    }

    if(!ancflg){
	cerr<<"Error: the root pop was not found"<<endl;
	return 1;
    }


    AlleleRecords * test;
    //unsigned int i;
    while(mp.hasData()){
	test = mp.getData();

        cout<<test->chr<<"\t"<<stringify(test->coordinate)<<"\t"<<test->ref<<","<<test->alt<<"\t";
	cout<<test->vectorAlleles->at(rootidx) <<"\t"<<test->vectorAlleles->at(ancidx) <<endl;
    }


    cerr<<"Program usePopAsRootAnc terminated gracefully"<<endl;

    return 0;
}

