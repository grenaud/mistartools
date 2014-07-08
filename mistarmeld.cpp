/*
 * mistarmeld.cpp
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

    bool keepOrig=false;

    if(argc == 1 ||
       (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
       ){
	cerr<<"usage: "<<argv[0]<<" <options> [mistar file zipped or normal] \"popToMerge1,popToMerge2,....\" \"newid\"\nwill print to stdout"
	    <<endl
	    <<"This program will merge different specified populations into a single one "<<endl
	    <<"ex: "<<argv[0]<<" data.mst \"Papuan,Austalian\" \"oceanians\""<<endl
	    <<"Options:"<<endl
	    <<"\t--keep\t\tKeep the original populations in the output (Default "+boolStringify(keepOrig)+" )\n"
	    // <<"\t--noprivate\t\tDo not allow private mutations (Default "+boolStringify(noprivate)+" )\n"

	    <<endl;
	
	return 1;
    }
    
    //all but two last
    for(int i=1;i<(argc-3);i++){ 
	

        if( string(argv[i]) == "--keep"  ){
	    keepOrig=true;
            continue;
        }

        // if( string(argv[i]) == "--noprivate"  ){
	//     noprivate=true;
        //     continue;
        // }

	cerr<<"Error unknown option "<<argv[i]<<endl;
	exit(1);
    }

    MistarParser mp                  (string(argv[argc-3]));
    string pop2mergeString           =string(argv[argc-2]);
    string mergedpopName             =string(argv[argc-1]);

    vector<string> tomerge = allTokens( pop2mergeString,',');
    set<unsigned int> indexPopTomerge;
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    cout<<"#MISTAR"<<endl;    	
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;
    cout<<"#MISTARMELD:"<<endl;
    //for(unsigned int i=0;i<vectorOfMP.size();i++){ 
    cout<<"#UNIONMELD#"<<endl;
    cout<<""<<mp.getHeader("#\t")<<"\n";
    //}

    cout<<"#chr\tcoord\tREF,ALT"<<"\t";


    for(unsigned int i=0;i<mp.getPopulationsNames()->size();i++){
	for(unsigned int j=0;j<tomerge.size();j++){
	    if( tomerge[j] == mp.getPopulationsNames()->at(i) ){
		indexPopTomerge.insert(i);
	    }
	}

    }

    if(indexPopTomerge.size() != tomerge.size()){
	cerr<<"Error: some of the populations to merge "<<pop2mergeString<<" were not found"<<endl;
	return 1;
    }

    // cout<<vectorToString(indexPopTomerge)<<endl;
    // return 1;

    // cout<<vectorToString( *(mp.getPopulationsNames())," ")<<endl;
    for(unsigned int i=0;i<mp.getPopulationsNames()->size();i++){
	if( indexPopTomerge.find(i) ==   indexPopTomerge.end() ){
	    cout<<mp.getPopulationsNames()->at(i)<<"\t";
	}else{
	    if(keepOrig)
		cout<<mp.getPopulationsNames()->at(i)<<"\t";
	    //skip
	}
    }

    cout<<mergedpopName<<endl;

    AlleleRecords * test;
    unsigned int i;
    while(mp.hasData()){
	// cout<<endl<<"data"<<endl;
	test = mp.getData();
	// cout<<test<<endl;
	cout<<test->chr<<"\t"<<stringify(test->coordinate)<<"\t"<<test->ref<<","<<test->alt<<"\t";

	SingleAllele newSingleAllele;
	// newSingleAllele.refCount=0;
	// newSingleAllele.altCount=0;
	// newSingleAllele.isCpg=false;

	for(i=0;i<test->vectorAlleles->size();i++){
	    if(indexPopTomerge.find(i) ==   indexPopTomerge.end() ){//print normally if not found
		cout<<test->vectorAlleles->at(i) <<"\t";
	    }else{
		if(keepOrig)
		    cout<<test->vectorAlleles->at(i) <<"\t";
		newSingleAllele+=test->vectorAlleles->at(i);
		// newSingleAllele.refCount += test->vectorAlleles->at(i).refCount;
		// newSingleAllele.altCount += test->vectorAlleles->at(i).altCount;
		// newSingleAllele.isCpg     = newSingleAllele.isCpg || test->vectorAlleles->at(i).isCpg; // if one is a CpG, the merged one is a CpG	       
	    }
	}
	cout<<newSingleAllele <<endl; 

    }


    cerr<<"Program mistarmeld terminated gracefully"<<endl;

    return 0;
}

