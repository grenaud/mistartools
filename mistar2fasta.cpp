/*
 * mistar2fasta
 * Date: Mar-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "MistarParser.h"

using namespace std;

int main (int argc, char *argv[]) {

    bool printRoot=true;
    bool produceTwo=false;

    string usage=string(""+string(argv[0])+" <options> [mistar file] "+
			"\n\nThis program takes a mistar matrix and prints a FASTA file using the allele information\nwith one record per population. Each site generates one base pair.\n\n"+
			"\tOptions\n"+			
			"\t\t"+"--noanc"+"\t"+"Do not print the root/anc (Default: "+boolStringify(printRoot)+" )\n"
			"\t\t"+"--het"+"\t"+"Produce two fasta files containing the alleles for het sites (Default: "+boolStringify(printRoot)+" )\n"		
	);
    //    cout<<argc<<endl;
    if(argc < 2 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    //all but last arg
    for(int i=0;i<(argc-1);i++){ 

	if( string(argv[i]) == "--noanc"){
	    printRoot=false;
	}

	if( string(argv[i]) == "--het"){
	    produceTwo=true;
	}

    }



    MistarParser mp         (argv[argc-1]);
        // cout<<"ok"<<endl;    

    AlleleRecords * record;
 
    unsigned int firstIndex=0;
    if(!printRoot)
	firstIndex=2;

    vector<string> deflines;
    vector<string> sequences;

    if(produceTwo){

	for(unsigned int i=firstIndex;i<mp.getPopulationsNames()->size();i++){
	    //indFileS<<mp.getPopulationsNames()->at(i)<<"\tU\t"<<mp.getPopulationsNames()->at(i)<<endl;	
	    deflines.push_back(mp.getPopulationsNames()->at(i)+"-1");
	    sequences.push_back("");
	    deflines.push_back(mp.getPopulationsNames()->at(i)+"-2");
	    sequences.push_back("");
	} 

    }else{
	for(unsigned int i=firstIndex;i<mp.getPopulationsNames()->size();i++){
	    //indFileS<<mp.getPopulationsNames()->at(i)<<"\tU\t"<<mp.getPopulationsNames()->at(i)<<endl;	
	    deflines.push_back(mp.getPopulationsNames()->at(i));
	    sequences.push_back("");
	} 
	// indFileS.close();
    }
    
    //unsigned int counter=0;
    while(mp.hasData()){
	record = mp.getData();
	// //	cout<<"ok"<<endl;    
	// cout<<*record<<endl;
	// // if(!isResolvedDNA(record->alt))
	// //     continue;

	// snpFileS<<"snp#"<<(counter++)<<"\t"<<record->chr<<"\t"<<stringify(double(record->coordinate)/double(1000000))<<"\t"<<stringify(record->coordinate)<<"\t"<<record->ref<<"\t"<<record->alt<<endl;
	
	unsigned int firstIndex=0;
	if(!printRoot)
	    firstIndex=2;
	unsigned int indexVec=0;
	if(produceTwo){

	    for(unsigned int i=firstIndex;i<record->vectorAlleles->size();i++){
		char c = record->vectorAlleles->at(i).generateRandomAllele(record->ref,record->alt);//otherwise, a file will have the major allele and the other one the minor.
		char c2;
		if(record->vectorAlleles->at(i).isHeterozygous()){
		    sequences[ int(indexVec*2.0)   ]  += c;

		    if(record->ref == c)
			c2=record->alt;
		    else
			c2=record->ref;
		    
		    sequences[ int(indexVec*2.0)+1 ]  += c2;

		}else{
		    sequences[ int(indexVec*2.0)   ]  += c;
		    sequences[ int(indexVec*2.0)+1 ]  += c;		
		}
		    
		indexVec++;
	    }
	    
	}else{
	    for(unsigned int i=firstIndex;i<record->vectorAlleles->size();i++){
		sequences[ indexVec++ ]  += record->vectorAlleles->at(i).generateRandomAlleleBias(record->ref,record->alt);
	    }
	}	
	// genoFileS<<endl;

    }
    // genoFileS.close();
    // snpFileS.close();

    for(unsigned int i=0;i<sequences.size();i++){
	cout<<">"<<deflines[i]<<endl;
	cout<<sequences[i]<<endl;

    }
    
    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}

