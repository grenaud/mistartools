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

    bool printRoot=true;

    string usage=string(""+string(argv[0])+" <options> [mistar file] [out (.bed)] [out (.bim)] [out SNP file (.fam)] "+
			"\n\nThis program takes a mistar matrix and prints the genotype and SNP file\n\n"+
			"\tOptions\n"+			
			"\t\t"+"--noanc"+"\t"+"Do not print the root/anc (Default: "+boolStringify(printRoot)+" )\n"
			);

    if(argc < 4 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    //all but last 4
    for(int i=0;i<(argc-4);i++){ 

	if( string(argv[i]) == "--noanc"){
	    printRoot=false;
	}

    }



    MistarParser mp         (argv[argc-4]);
    string bedFile  = string(argv[argc-3]);
    string bimFile  = string(argv[argc-2]);
    string famFile  = string(argv[argc-1]);

    ofstream bedFileS;
    ofstream bimFileS; 
    ofstream famFileS; 
	
    bedFileS.open(bedFile.c_str(),   ios::out| ios::binary);
    bimFileS.open(bimFile.c_str(),   ios::out);
    famFileS.open(famFile.c_str(),   ios::out);

    
    if (!bedFileS.good()){  cerr << "Unable to open file "<<bedFile<<endl; return 1;   }
    if (!bimFileS.good()){  cerr << "Unable to open file "<<bimFile<<endl; return 1;   }
    if (!famFileS.good()){  cerr << "Unable to open file "<<famFile<<endl; return 1;   }
    

    AlleleRecords * record;
 
    unsigned int firstIndex=0;

    if(!printRoot)
	firstIndex=2;

    for(unsigned int i=firstIndex;i<mp.getPopulationsNames()->size();i++){
    	famFileS<<mp.getPopulationsNames()->at(i)<<"\t"<<mp.getPopulationsNames()->at(i)<<"\t"<<mp.getPopulationsNames()->at(i)<<"\t"<<mp.getPopulationsNames()->at(i)<<"\t1\t-9"<<endl;	
    } 
    famFileS.close();

    //magic number
    char c = 108;
    bedFileS.write( (char *)&c, sizeof(c));
    c = 27;
    bedFileS.write( (char *)&c, sizeof(c));

    //snp major
    c = 1;
    bedFileS.write( (char *)&c, sizeof(c));



    unsigned int counter=0;
    while(mp.hasData()){
	record = mp.getData();
	
	if(!isResolvedDNA(record->alt))//skip non-seg sites
	    continue;

	unsigned int firstIndex=0;
	if(!printRoot)
	    firstIndex=2;


	if(!printRoot){
	    unsigned int refSum=0;
	    unsigned int altSum=0;
	    
	    
	    for(unsigned int i=firstIndex;i<record->vectorAlleles->size();i++){
		refSum+=record->vectorAlleles->at(i).getRefCount() ;
		altSum+=record->vectorAlleles->at(i).getAltCount() ;		
	    }

	    if( (altSum==0) || (refSum==0) )//skip non-seg sites in all but the root/anc
		continue;

	    
	}
	    
	
	bimFileS<<record->chr<<"\t"<<"snp#"<<(counter++)<<"\t"<<stringify(double(record->coordinate)/double(1000000))<<"\t"<<stringify(record->coordinate)<<"\t"<<record->ref<<"\t"<<record->alt<<endl;
	
	

	char byteToWrite=0; //data to write
	int storedInByte=0; //how much data is there in byteToWrite
	
	for(unsigned int i=firstIndex;i<record->vectorAlleles->size();i++){	    

	    char   toStore = record->vectorAlleles->at(i).printBinaryPLINK();
	    //cout<<i<<" "<<var2binary(toStore)<<endl;
	    if( (storedInByte%4) == 3){//to write
		byteToWrite |= (  toStore<< ((storedInByte%4)*2) );		
		//cout<<"snp# "<<(counter-1)<<" "<<var2binary(byteToWrite)<<endl;
		bedFileS.write( (char *)&byteToWrite, sizeof(byteToWrite));		
		byteToWrite  = 0;
		storedInByte = 0;
	    }else{
		byteToWrite |= (  toStore<< ((storedInByte%4)*2) );
		storedInByte++;
	    }
	    
	}

	if(storedInByte!=0){
	    //cout<<"padpre snp "<<(counter-1)<<" "<<var2binary(byteToWrite)<<endl;  
	    for(int k=3;k>=storedInByte;k--){
		char zeroC= 3; //missing data
		byteToWrite |= (  zeroC << ((k%4)*2) );		
	    }
	    //cout<<"padpos snp "<<(counter-1)<<" "<<var2binary(byteToWrite)<<endl;  
	    bedFileS.write( (char *)&byteToWrite, sizeof(byteToWrite));		
	}
	
	//bedFileS<<endl;

    }
    bedFileS.close();
    bimFileS.close();
    
    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}

