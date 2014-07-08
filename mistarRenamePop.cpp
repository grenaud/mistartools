/*
 * mistarrenamepop.cpp
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
    // cerr<<"arg "<<endl;
    
    if(argc == 1 ||
       (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
       ){
	cerr<<"usage: "<<argv[0]<<" <options> [mistar file] \"popOldName1,popOldName2,....\" \"popNewName1,popNewName2,....\"\nwill print to stdout"
	    <<endl
	    <<"This program will rename different specified populations  "<<endl
	    <<"ex: "<<argv[0]<<" data.mst \"Papuan,Austalian\" \"Oceanians1,Oceanians2\""<<endl
	    <<"Options:"<<endl
	    // <<"\t--justtransv\t\tOnly allow transversions (Default "+boolStringify(limitToTransversions)+" )\n"
	    // <<"\t--noprivate\t\tDo not allow private mutations (Default "+boolStringify(noprivate)+" )\n"

	    <<endl;	
	return 1;
    }
    
    // cerr<<"arg "<<endl;

    
    // //all but two last
    // for(int i=1;i<(argc);i++){ 
	
    // 	cerr<<"arg "<<argc<<"\t"<<i<<"\t"<<string(argv[i])<<endl;
    // }
    //     // if( string(argv[i]) == "--justtransv"  ){
    // 	//     limitToTransversions=true;
    //     //     continue;
    //     // }

    //     // if( string(argv[i]) == "--noprivate"  ){
    // 	//     noprivate=true;
    //     //     continue;
    //     // }

    // 	cerr<<"Error unknown option "<<argv[i]<<endl;
    // 	exit(1);
    // }

    MistarParser mp                  (string(argv[argc-3]));
    string pop2oldString           =string(argv[argc-2]);
    string pop2newString           =string(argv[argc-1]);
    set<unsigned int> indexPopToRename;
    cerr<<"replacing "<<pop2oldString<<" to "<<pop2newString<<" in file "<<string(argv[argc-3])<<endl;
    vector<string> popold = allTokens( pop2oldString,',');
    vector<string> popnew = allTokens( pop2newString,',');

    if(popnew.size() != popold.size()){
	cerr<<"Error: the number of populations from the old names and new does not jibe"<<endl;
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
    cout<<"#MISTARRENAMEPOP:"<<endl;
    //for(unsigned int i=0;i<vectorOfMP.size();i++){ 
    cout<<"#RENAMEPOP#"<<endl;
    cout<<""<<mp.getHeader("#\t")<<"\n";
    //}

    cout<<"#chr\tcoord\tREF,ALT"<<"\t";
    vector<string> newPopVector;

    for(unsigned int i=0;i<mp.getPopulationsNames()->size();i++){
	//cerr<<mp.getPopulationsNames()->at(i)<<endl;
	bool found=false;
	for(unsigned int j=0;j<popold.size();j++){
	    if( popold[j] == mp.getPopulationsNames()->at(i) ){
		if(found){
		    cerr<<"Error: some of the populations  "<<popold[j]  <<" were found twice"<<endl;
		    return 1;
		}
		
		newPopVector.push_back( popnew[j]  );
		indexPopToRename.insert(i);

		found=true;
	    }
	}

	if(!found){
	    newPopVector.push_back(  mp.getPopulationsNames()->at(i) );
	}
    }

    if(indexPopToRename.size() != popold.size()){
	cerr<<"Error: some of the populations to merge "<< pop2oldString <<" were not found"<<endl;
	return 1;
    }

    cout<<vectorToString(newPopVector,"\t")<<endl;
    // return 1;

    // cout<<vectorToString( *(mp.getPopulationsNames())," ")<<endl;

    AlleleRecords * test;
    //unsigned int i;
    while(mp.hasData()){
	 // cout<<endl<<"data"<<endl;
	test = mp.getData();
	cout<<*test<<endl;
	// cout<<test->chr<<"\t"<<stringify(test->coordinate)<<"\t"<<test->ref<<","<<test->alt<<"\t";

	// SingleAllele newSingleAllele;
	// // newSingleAllele.refCount=0;
	// // newSingleAllele.altCount=0;
	// // newSingleAllele.isCpg=false;

	// for(i=0;i<test->vectorAlleles->size();i++){
	//     if(indexPopTomerge.find(i) ==   indexPopTomerge.end() ){//print normally if not found
	// 	cout<<test->vectorAlleles->at(i) <<"\t";
	//     }else{
	// 	newSingleAllele+=test->vectorAlleles->at(i);
	// 	// newSingleAllele.refCount += test->vectorAlleles->at(i).refCount;
	// 	// newSingleAllele.altCount += test->vectorAlleles->at(i).altCount;
	// 	// newSingleAllele.isCpg     = newSingleAllele.isCpg || test->vectorAlleles->at(i).isCpg; // if one is a CpG, the merged one is a CpG	       
	//     }
	// }
	// cout<<newSingleAllele <<endl; 

    }


    cerr<<"Program mistarmeld terminated gracefully"<<endl;

    return 0;
}

