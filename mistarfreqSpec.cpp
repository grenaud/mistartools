/*
 * mistarfreqspec
 * Date: May-28-2013 
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
    bool useAnc =false;
    bool useRoot=false;

    string usage=string(""+string(argv[0])+" <options> [mistar file]"+
			"\nThis program will print the number of observed alleles for the reference and alternative alleles:\n\n"+
			"\t\t"+"Options:\n"+
			"\t\t"+"--useanc\t\tUse the ancestral allele to report the frequency\n"+
			"\t\t"+"--useroot\t\tUse the root allele to report the frequency\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    
    
    //starts at 1 and except the last two
    for(int i=1;i<(argc-1);i++){ 
	if(string(argv[i]) == "--useanc" ) {
	    useAnc  = true;
	    continue;
	}

	if(string(argv[i]) == "--useroot" ) {
	    useRoot = true;
	    continue;
	}


        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1;           

    }


    if(useRoot && useAnc){
	cerr<<"Cannot use both the ancestor and the root"<<endl;
        return 1;           
    }

    unsigned int totalRecords=0;

    MistarParser mp (argv[argc-1]);
    AlleleRecords * dataRow;


    
    while(mp.hasData()){
	dataRow = mp.getData();
	totalRecords++;

	//non-seg site, no point in looking at those
	if(!isResolvedDNA(dataRow->alt))
	    continue;

	unsigned int refCounter=0;
	unsigned int altCounter=0;
	bool rootIsRef = false;
	bool ancIsRef  = false;

	if(useRoot){
	    if(dataRow->vectorAlleles->at(0).getRefCount() == 0 &&
	       dataRow->vectorAlleles->at(0).getAltCount() == 0 ){
		continue;	    
	    }
	    

	    if(dataRow->vectorAlleles->at(0).getAltCount() != 0 ){
		if(dataRow->vectorAlleles->at(0).getRefCount() != 0 ){
		    cerr<<"Cannot determine the root allele for "<<*dataRow<<endl;
		    return 1;      
		}	
		rootIsRef=false;
	    }else{
		rootIsRef=true;
	    }
    
	}

	if(useAnc){
	    if(dataRow->vectorAlleles->at(1).getRefCount() == 0 &&
	       dataRow->vectorAlleles->at(1).getAltCount() == 0 ){
		continue;
	    }
	    
	    if(dataRow->vectorAlleles->at(1).getAltCount() != 0 ){
		if(dataRow->vectorAlleles->at(1).getRefCount() != 0 ){
		    cerr<<"Cannot determine the root allele for "<<*dataRow<<endl;
		    return 1;      
		}	
		ancIsRef=false;
	    }else{
		ancIsRef=true;
	    }
	}




	for(unsigned j=2;j<dataRow->vectorAlleles->size();j++){
	    //undefined site
	    refCounter+=dataRow->vectorAlleles->at(j).getRefCount();
	    altCounter+=dataRow->vectorAlleles->at(j).getAltCount();	
	}


	cout<<dataRow->coordinate<<"\t";
	if(useRoot){
	    if(rootIsRef)
		cout<<refCounter<<"\t"<<altCounter<<endl;
	    else
		cout<<altCounter<<"\t"<<refCounter<<endl;
	    continue;
	}

	if(useAnc){
	    if(ancIsRef)
		cout<<refCounter<<"\t"<<altCounter<<endl;
	    else
		cout<<altCounter<<"\t"<<refCounter<<endl;
	    continue;
	}

	cout<<refCounter<<"\t"<<altCounter<<endl;

    }

    cerr<<"Program "<<argv[0]<<" looked at  "<<totalRecords<<"records, terminated gracefully";

    return 0;
}

