/*
 * mistar2AlleleMatrix
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

    bool printRoot=false;

    string usage=string(""+string(argv[0])+" <options> [mistar file] [out genotype file (.geno)] [out SNP file (.snp)] [out SNP file (.ind)] "+
			"\n\nThis program takes a mistar matrix and prints the genotype and SNP file\n\n"+
			"\tOptions\n"+			
			"\t\t"+"--withanc"+"\t"+"Print the root/anc (Default: "+boolStringify(printRoot)+" )\n"
			);

    if(argc < 4 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    //all but last 4
    for(int i=0;i<(argc-4);i++){ 

	if( string(argv[i]) == "--withanc"){
	    printRoot=true;
	}

    }



    MistarParser mp         (argv[argc-4]);
    string genoFile = string(argv[argc-3]);
    string snpFile  = string(argv[argc-2]);
    string indFile  = string(argv[argc-1]);

    ofstream genoFileS;
    ofstream snpFileS; 
    ofstream indFileS; 
	
    genoFileS.open(genoFile.c_str(), ios::out);
    snpFileS.open(snpFile.c_str(),   ios::out);
    indFileS.open(indFile.c_str(),   ios::out);

    
    if (!genoFileS.good()){       cerr << "Unable to open file "<<genoFile<<endl;       return 1;     }
    if (!snpFileS.good()){        cerr << "Unable to open file "<<snpFile<<endl;        return 1;     }
    if (!indFileS.good()){        cerr << "Unable to open file "<<indFile<<endl;        return 1;     }
    

    AlleleRecords * record;
 
    unsigned int firstIndex=0;
    if(!printRoot)
	firstIndex=2;
    for(unsigned int i=firstIndex;i<mp.getPopulationsNames()->size();i++){
	indFileS<<mp.getPopulationsNames()->at(i)<<"\tU\t"<<mp.getPopulationsNames()->at(i)<<endl;
	
    } 
    indFileS.close();

    
    unsigned int counter=0;
    while(mp.hasData()){
	record = mp.getData();
	
	if(!isResolvedDNA(record->alt))
	    continue;

	if( (counter%1000000) ==0 ){
	    cerr<<"Produced "<<counter<<" sites at "<<record->chr<<":"<<record->coordinate<<endl;
	}
	snpFileS<<"snp#"<<(counter++)<<"\t"<<record->chr<<"\t"<<stringify(double(record->coordinate)/double(1000000))<<"\t"<<stringify(record->coordinate)<<"\t"<<record->ref<<"\t"<<record->alt<<endl;
	
	unsigned int firstIndex=0;
	if(!printRoot)
	    firstIndex=2;

	for(unsigned int i=firstIndex;i<record->vectorAlleles->size();i++){
	    genoFileS<<record->vectorAlleles->at(i).printEIGENSTRAT();	   
	} 
	genoFileS<<endl;

    }
    genoFileS.close();
    snpFileS.close();
    
    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}

