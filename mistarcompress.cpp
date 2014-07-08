/*
 * mistarcompress
 * Date: Apr-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "MistarParser.h"
#include "mistarOperations.h"

using namespace std;


int main (int argc, char *argv[]) {
    //    bool printChr=false;
    string fastaIndex  =  "/mnt/454/Altaiensis/users/gabriel/faidx/index.hg19.fai"  ;

    string usage=string(""+string(argv[0])+" <options> [mistar file] "+
			"\n\nThis program takes a mistar matrix (zipped or not) and prints it as binary\n\n"+
			"\tOptions\n"+
			// "\t\t"+"--chrcoord"+"\t"+"Print the chr/coord (Default: "+stringify(printChr)+" )"+
			"\t\t"+"--fai     [file]"            + "\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: "+fastaIndex+")\n"+

			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }



    //starts at 1 and except the last two
    for(int i=1;i<(argc-2);i++){ 

	if(string(argv[i]) == "--fai"){
	    fastaIndex=string(argv[i+1]);
	    i++;
	    continue;
	}

        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1;           

    }

    unsigned int totalRecords=0;

    vector<chrinfo> chrFound;
    uint64_t genomeLength;
    readFastaIndex(fastaIndex,chrFound,genomeLength);

    MistarParser mp       (argv[argc-1]);

    char magicstr [4] = {'M','S','T','\1'};
    
    //magicstr="MST";
    if(write(1,&magicstr,sizeof(magicstr)) == -1 ){   cerr<<"Write error"<<endl; 	    return 1; 	}     

    

    stringstream toflush;
    //cout<<""<<mp.getHeaderNoDefline()<<"\n";
    toflush<<mp.getHeaderNoDefline()<<endl;

    map<string,uint16_t> chr2index;
    uint16_t     chrCurrentIndex=0;

    for(unsigned j=0;j<(chrFound.size());j++){
	//cout<<"#SQ\t"<<"SN:"<<chrFound[j].name<<"\tLN:"<<chrFound[j].length<<endl;
	toflush<<"#SQ\t"<<"SN:"<<chrFound[j].name<<"\tLN:"<<chrFound[j].length<<endl;
	//toflush+=mp.getHeaderNoDefline();

	chr2index[chrFound[j].name]=chrCurrentIndex++;
	if(chrCurrentIndex == 0xFFFF){
	    cerr<<"Too many chromosomes for this build, more than 65535"<<endl;
	    return 1;
	}
    }
    
    toflush<<mp.getDefline()<<endl;


    
    uint32_t sizeHeader= toflush.str().size();
    

    if(write(1,&sizeHeader,sizeof(sizeHeader)) == -1 ){   cerr<<"Write error"<<endl; 	    return 1; 	} 
    char  toflush_cptr [sizeHeader];
    strcpy(toflush_cptr,toflush.str().c_str());
    //cerrr<<
    for(uint32_t i=0;i<sizeHeader;i++){
        char towrite=char(toflush_cptr[i]);
	if(write(1,&towrite,sizeof(towrite)) == -1 ){   cerr<<"Write error"<<endl; 	    return 1; 	} 
    }
    

    uint32_t sizePops=mp.getPopulationsNames()->size();

    if(write(1,&sizePops,sizeof(sizePops)) == -1 ){   cerr<<"Write error"<<endl; 	    return 1; 	} 
    
    
    //compute size of data row
    //chr uin16_t 2
    //coordinate uin32_t 4
    //char ref 1
    //char alt 1
    //short,short,char X  mp.getPopulationsNames()->size()  5
    //total:
    // 2+4+1+1+ 5xmp.getPopulationsNames()->size() 
    // 8 + 5 x sizePop
    
    

    AlleleRecords * dataRow;
    char  tempCh;
    short tempSh;
    while(mp.hasData()){
	dataRow = mp.getData();
	if(chr2index.find(dataRow->chr) == chr2index.end()){
	    cerr<<"Cannot find chr "<<dataRow->chr<<" in index "<<endl;
	    return 1;
	}
	
	//write index of chr 2 bytes
	if(write(1,&chr2index[dataRow->chr],sizeof(chr2index[dataRow->chr])) == -1 ) {   cerr<<"Write error"<<endl; 	    return 1; 	}
	//write coordinate 4 bytes
	if(write(1,&dataRow->coordinate,    sizeof(dataRow->coordinate))  == -1 )    {   cerr<<"Write error"<<endl; 	    return 1; 	}

	//ref 1 byte
	tempCh= char(base2int(dataRow->ref));
	if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl; 	    return 1; 	}
	//alt 1 byte
	tempCh= char(base2int(dataRow->alt));
	if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl; 	    return 1; 	}
	

	for(unsigned j=0;j<(mp.getPopulationsNames()->size());j++){
	    //2 bytes
	    tempSh= short( dataRow->vectorAlleles->at(j).getRefCount() );
	    if(write(1,&tempSh,sizeof(tempSh))  == -1 ){   cerr<<"Write error"<<endl; 	    return 1; 	}
	    //2 bytes
	    tempSh= short( dataRow->vectorAlleles->at(j).getAltCount() );
	    if(write(1,&tempSh,sizeof(tempSh))  == -1 ){   cerr<<"Write error"<<endl; 	    return 1; 	}
	    //1 byte
	    tempCh= short( dataRow->vectorAlleles->at(j).getIsCpg() );
	    if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl; 	    return 1; 	}
	}
	//todelete 
	// uint16_t     placeholder= 0xFFFF;
	// if(write(1,&placeholder,sizeof(uint16_t))  == -1 ){   cerr<<"Write error"<<endl;    return 1; 	}

	totalRecords++;
    }

    cerr<<"Program "<<argv[0]<<" wrote "<<totalRecords<<" terminated gracefully";

    return 0;
}

