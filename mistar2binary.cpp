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
    bool printChr=false;
    string usage=string(""+string(argv[0])+" <options> [mistar file] [comma separated group] "+
			"\n\nThis program takes a mistar matrix and prints the alleles\nas a binary matrix (0=ancestral,1=derived)\n\n"+
			"\tOptions\n"+
			"\t\t"+"--chrcoord"+"\t"+"Print the chr/coord (Default: "+stringify(printChr)+" )\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }



    //starts at 1 and except the last two
    for(int i=1;i<(argc-2);i++){ 
	if(string(argv[i]) == "--chrcoord" ) {
	    printChr=true;
	    continue;
	}

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
    string g1 = string(argv[argc-1]);
    //string g2 = string(argv[argc-1]);
    vector<string> g1v = allTokens(g1,',');
    //vector<string> g2v = allTokens(g2,',');
    vector<unsigned int>   g1i;
    vector<string>   g1s;

    //vector<unsigned int>   g2i;

    if(printChr){
	cout<<"#chr\tcoord\t";
    }else{
	cout<<"#";
    }
	    

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
    cout<<vectorToString(g1s,"\t")<<endl;
    // exit(1);

    AlleleRecords * dataRow;
    // cout<<"#MISTARFILTER:nosharing "<<g1<<"-"<<g2<<endl;
    // cout<<""<<mp.getHeader("#\t")<<"\n";
    // cout<<""<<mp.getDefline()<<"\n";

    while(mp.hasData()){
	dataRow = mp.getData();
	totalRecords++;
	vector<bool> alleleBinToPrint;
	char ancestralAllele;
	//getting ancestral allele

	if(dataRow->vectorAlleles->at(1).getRefCount() == 0 && 
	   dataRow->vectorAlleles->at(1).getAltCount() == 0 )
	    goto nextiteration;
	
	ancestralAllele=sampleRandomRefAltAllele(dataRow->ref,
						 dataRow->alt,
						 dataRow->vectorAlleles->at(1).getRefCount(), 
						 dataRow->vectorAlleles->at(1).getAltCount());

	
	for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){
	    
	    
	    for(unsigned k=0;k<g1i.size();k++){
		if(j==g1i[k]){
		    //skip undefined sites 
		    if(dataRow->vectorAlleles->at(j).getRefCount() == 0 && 
		       dataRow->vectorAlleles->at(j).getAltCount() == 0 )
			goto nextiteration;			       
		    
		    if(ancestralAllele==sampleRandomRefAltAllele(dataRow->ref,
								 dataRow->alt,
								 dataRow->vectorAlleles->at(j).getRefCount(), 
								 dataRow->vectorAlleles->at(j).getAltCount())){
			alleleBinToPrint.push_back(false);
		    }else{
			alleleBinToPrint.push_back(true);
		    }

		}
	    }
		

		    

	}//end for each column of data
		

	if(printChr){
	    cout<<dataRow->chr<<"\t"<<dataRow->coordinate<<"\t";
	}

	cout<<vectorToString(alleleBinToPrint,"\t")<<endl;
	keptRecords++;
    nextiteration:	    
	continue;
    }

    cerr<<"Program "<<argv[0]<<" wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully";

    return 0;
}

