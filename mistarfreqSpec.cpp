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
    bool onlysegsite =false;
    bool splitpop =false;
    bool usefreq =false;

    string usage=string(""+string(argv[0])+" <options> [mistar file]"+
			"\nThis program will print the number of observed alleles for the reference and alternative alleles\nThe left column is the reference count. The order can be changed using --useanc and --useroot options:\n\n"+
			"\t\t"+"Options:\n"+
			"\t\t"+"--splitpop\t\tSplit pop.\n"+
			"\t\t"+"--freq\t\tOutput frequencies\n"+

			"\t\t"+"--onlysegsite\t\tUse only segregating sites\n"+
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
	if(string(argv[i]) == "--onlysegsite" ) {
	    onlysegsite = true;
	    continue;
	}

	if(string(argv[i]) == "--freq" ) {
	    usefreq= true;
	    continue;
	}

	if(string(argv[i]) == "--splitpop" ) {
	    splitpop  = true;
	    continue;
	}

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
	if(onlysegsite)
	    if(!isResolvedDNA(dataRow->alt))
		continue;

	unsigned int refCounter=0;
	unsigned int altCounter=0;
	vector< pair<unsigned int,unsigned int> > alleleC;

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



	if(splitpop){
	    for(unsigned j=2;j<dataRow->vectorAlleles->size();j++){
		alleleC.push_back( make_pair<unsigned int,unsigned int>( dataRow->vectorAlleles->at(j).getRefCount(),
									 dataRow->vectorAlleles->at(j).getAltCount() ) );
	    }

	    cout<<dataRow->coordinate;
	    bool flipAlt=false;
	    if(useRoot){
		if(!rootIsRef)
		    flipAlt=true;
	    }
	    
	     if(useAnc){
		if(!ancIsRef)
		    flipAlt=true;
	     }
	     
	     if(flipAlt){
		 
		 for(unsigned j=0;j<alleleC.size();j++){
		     double total=1.0;
		     if(usefreq)
			 total= double(alleleC[j].first) + double(alleleC[j].second);
		     cout<<"\t"<<alleleC[j].second/total<<"\t"<<alleleC[j].first/total;

		 }
		 
	     }else{
		 for(unsigned j=0;j<alleleC.size();j++){
		     double total=1.0;
		     if(usefreq)
			 total= double(alleleC[j].first) + double(alleleC[j].second);
		     cout<<"\t"<<alleleC[j].first/total<<"\t"<<alleleC[j].second/total;
		     
		 }
	     }
	     cout<<endl;
	}else{
	    for(unsigned j=2;j<dataRow->vectorAlleles->size();j++){
		//undefined site
		refCounter+=dataRow->vectorAlleles->at(j).getRefCount();
		altCounter+=dataRow->vectorAlleles->at(j).getAltCount();	
	    }

	    cout<<dataRow->chr<<"\t"<<dataRow->coordinate<<"\t";
	    double total=1.0;
	    if(usefreq){
		total= double(refCounter) + double(altCounter);
	    }

	    if(useRoot){
		if(rootIsRef)
		    cout<<double(refCounter)/total<<"\t"<<double(altCounter)/total<<endl;
		else
		    cout<<double(altCounter)/total<<"\t"<<double(refCounter)/total<<endl;
		continue;
	    }

	    if(useAnc){
		if(ancIsRef)
		    cout<<double(refCounter)/total<<"\t"<<double(altCounter)/total<<endl;
		else
		    cout<<double(altCounter)/total<<"\t"<<double(refCounter)/total<<endl;
		continue;
	    }

	    cout<<double(refCounter)/total<<"\t"<<double(altCounter)/total<<endl;
	}



    }

    cerr<<"Program "<<argv[0]<<" looked at  "<<totalRecords<<"records, terminated gracefully";

    return 0;
}

