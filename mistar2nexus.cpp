/*
 * mistar2next
 * Date: May-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "MistarParser.h"

using namespace std;

int main (int argc, char *argv[]) {
    bool printOnlySeq =false;
    bool allDefined=false;
    
    string usage=string(""+string(argv[0])+" <options> [mistar file] "+
			"\n\nThis program takes a mistar matrix and prints the alleles\nin nexus format\n\n"+
			"\tOptions:\n"+
			"\t\t"+"--group"+"\t"+"Just use these populations, comma separated list (Default: use all )\n"+
			"\t\t"+"--nonseg"+"\t"+"Print only segregating sites (Default: "+booleanAsString(printOnlySeq)+" )\n"
			"\t\t"+"--alldef"+"\t"+"Print only sites that are defined in the groups (Default: "+booleanAsString(allDefined)+" )\n"
			);

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    string g1 = "";

    //starts at 1 and except the last two
    for(int i=1;i<(argc-2);i++){ 
	if(string(argv[i]) == "--group" ) {
	    g1=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "--nonseg" ) {
	    printOnlySeq=true;
	    continue;
	}

	if(string(argv[i]) == "--alldef" ) {
	    allDefined=true;
	    continue;
	}



        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1;           

    }

    unsigned int totalRecords=0;
    unsigned int keptRecords=0;

    // string programLine;
    // for(int i=0;i<(argc);i++){ 
    // 	programLine+=(string(argv[i])+" ");
    // }

    // cout<<"#MISTAR"<<endl;    	
    // cout<<"#PG:"<<programLine<<endl;
    // cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    // cout<<"#DATE: "<<getDateString()<<endl;

    MistarParser mp   (argv[argc-1]);
    //string g1 = string(argv[argc-1]);
    //string g2 = string(argv[argc-1]);

    if(g1==""){//empty
	g1=vectorToString( *(mp.getPopulationsNames()),",");
    }
    // cout<<g1<<endl;
    // return 1;

    vector<string> g1v = allTokens(g1,',');
    //vector<string> g2v = allTokens(g2,',');
    vector<unsigned int>   g1i;
    vector<string>   g1s;

    //vector<unsigned int>   g2i;

    // if(printChr){
    // 	cout<<"#chr\tcoord\t";
    // }else{
    // 	cout<<"#";
    // }
	    

    for(unsigned k=0;k<g1v.size();k++){
	bool found=false;
	for(unsigned i=0;i<(mp.getPopulationsNames()->size());i++){
	    if(mp.getPopulationsNames()->at(i) == g1v[k]){
		g1i.push_back(i);
		//g1s.push_back( mp.getPopulationsNames()->at(i) );
		found=true;
		break;
	    } 
	}
	if(!found){
	    cerr<<"Cannot find population "<<g1v[k]<<endl;
	    return 1;
	}
    }
    
    for(unsigned j=0;j<(mp.getPopulationsNames()->size());j++){
	for(unsigned k=0;k<g1i.size();k++){
	    if(j==g1i[k]){
		g1s.push_back( mp.getPopulationsNames()->at(j) );
	    }
	}
    }


    // cout<<vectorToString(g1i)<<endl;
    //cout<<vectorToString(g1s,"\t")<<endl;
    // exit(1);
    cout<<"#NEXUS"<<endl;
    cout<<"begin data;"<<endl;
    AlleleRecords * dataRow;
    vector<string> alllelesToPrint (g1s.size(),"") ;


    
    while(mp.hasData()){
	dataRow = mp.getData();
	totalRecords++;

	bool skipThis=false;
	bool foundRef=false;
	bool foundAlt=false;
	//doing basic filtering
	
	for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){

	    for(unsigned k=0;k<g1i.size();k++){
		if(j==g1i[k]){

		    //put ? for undefined sites
		    if(allDefined){
			if(dataRow->vectorAlleles->at(j).getRefCount() == 0 && 
			   dataRow->vectorAlleles->at(j).getAltCount() == 0 ){
			    skipThis=true;
			    break;
			}
		    }
		    
		    if(printOnlySeq){
			foundRef = foundRef || (dataRow->vectorAlleles->at(j).getRefCount() != 0);
			foundAlt = foundAlt || (dataRow->vectorAlleles->at(j).getAltCount() != 0);
		    }
		    
					   		    
		}
	    }

	    if(skipThis)
		break;
	}


	if(printOnlySeq){
	    if(!(foundRef && foundAlt))
		continue;
	}

	if(skipThis)
	    continue;
	// cout<<"foundRef "<<foundRef<<endl;
	// cout<<"foundAlt "<<foundAlt<<endl;

	//cout<<dataRow->coordinate<<endl;
	// return 1;
	for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){

	    for(unsigned k=0;k<g1i.size();k++){
		if(j==g1i[k]){

		    //put ? for undefined sites
		    if(dataRow->vectorAlleles->at(j).getRefCount() == 0 && 
		       dataRow->vectorAlleles->at(j).getAltCount() == 0 ){
			alllelesToPrint[k]+="?";
			continue;
		    }
			
		    
		    //goto nextiteration;			       
		    alllelesToPrint[k]+=sampleRandomRefAltAllele(dataRow->ref,
								 dataRow->alt,
								 dataRow->vectorAlleles->at(j).getRefCount(), 
								 dataRow->vectorAlleles->at(j).getAltCount());

		    
		}
	    }
		
	}
	keptRecords++;	

	    

    }//end for each column of data
		

    

    cout<<"dimensions ntax="<<g1s.size()<<" nchar="<<keptRecords<<";"<<endl;
    cout<<"format datatype=dna interleave=no missing=? gap=-;"<<endl;
    cout<<"matrix"<<endl<<endl;
    unsigned int maxnamesize=g1s[0].size();
    for(unsigned k=1;k<g1i.size();k++){
	if(g1s[k].size() > maxnamesize)
	    maxnamesize = g1s[k].size();
    }

    for(unsigned k=0;k<g1i.size();k++){
	cout<<g1s[k] <<string(maxnamesize-g1s[k].size(),' ')<<"   "<<alllelesToPrint[k]<<endl;
    }
    cout<<endl<<";"<<endl;
    cout<<"end;"<<endl;

    cerr<<"Program "<<argv[0]<<" wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

