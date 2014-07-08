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
			"\n\nThis program takes a binary compressed mistar matrix and prints it as a text file\n\n"+
			"\tOptions\n"+
			// "\t\t"+"--chrcoord"+"\t"+"Print the chr/coord (Default: "+stringify(printChr)+" )"+
			//"\t\t"+"--fai     [file]"            + "\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: "+fastaIndex+")\n"+

			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }



    //starts at 1 and except the last two
    for(int i=1;i<(argc-2);i++){ 

	// if(string(argv[i]) == "--fai"){
	//     fastaIndex=string(argv[i+1]);
	//     i++;
	//     continue;
	// }

        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1;           

    }

    //myFilezipped=new igzstream();
    igzstream myFilezipped ;
    string filename = string( argv[argc-1] );
    
    myFilezipped.open(filename.c_str(), ios::in);

   if (!myFilezipped.good()){
       cerr <<"Unable to open file : "<<filename<<endl;
       return 1;
    }

    const string magicNumber="MST\1";

    unsigned char  stringHeader [4];
    const unsigned int     lengthHeader  = 4;

    myFilezipped.read((char*)&stringHeader,lengthHeader);

    for(unsigned int i=0;i<lengthHeader;i++){
        if(magicNumber[i] != stringHeader[i]){
            cerr<<"MISTAR file does not begin with "<<magicNumber<<endl;
            exit(1);
        }
    }

    cout<<"#MISTAR"<<endl;      
    


    uint32_t sizeHeader;
    myFilezipped.read((char*)&sizeHeader,sizeof(sizeHeader));
    string toPrintHeader="";

    for(uint32_t i=0;i<sizeHeader;i++){
        char toread;
        myFilezipped.read((char*)&toread,sizeof(char));

	toPrintHeader+=toread;
    }

    unsigned int totalRecords=0;

    istringstream f (toPrintHeader);
    string line;    
    vector<string> chrKnown;
    
    while (getline(f, line)) {
	if(strBeginsWith(line,"#SQ")){
	    vector<string> tokensf = allTokens(line,'\t');
	    chrKnown.push_back(tokensf[1].substr(3));
	}else{
	    cout << line << std::endl;
	}
    }
    //cout<<vectorToString(chrKnown)<<endl;

    //return 1;
    uint32_t sizePops;
    myFilezipped.read((char*)&sizePops,sizeof(sizePops));
    //cout<<sizePops<<endl;

    size_t sizeRecord = 8+ 5*sizePops;
    char buffer [sizeRecord];
    while(!myFilezipped.eof()){
	myFilezipped.read((char*)&buffer,        sizeof(buffer));

	if(!myFilezipped){
	    //I have no clue why eof does not work
	    if(myFilezipped.gcount() != 0){
		cerr << "ERROR in decompression at record#"<< totalRecords<<", only: " << myFilezipped.gcount() << " could be read";
		return 1;
	    }else{
		break;
	    }

	}

	uint16_t chr; //2
	memcpy((char*)&chr,           buffer+0,    sizeof(chr));
	uint32_t coordinate ; //4
	memcpy((char*)&coordinate,    buffer+2,    sizeof(coordinate));
	char ref; //1
	memcpy((char*)&ref,           buffer+6,    sizeof(ref));
	char alt; //1
	memcpy((char*)&alt,           buffer+7,    sizeof(alt));

	cout<<chrKnown[chr]<<"\t";
	cout<<coordinate<<"\t";
	cout<<"NACGT"[ref]<<","<<"NACGT"[alt]<<"\t";

	for(unsigned j=0;j<(sizePops);j++){
	    short refC; //2
	    short altC; //2
	    char  cpgC; //1

	    memcpy((char*)&refC,       buffer+8 +5*j, sizeof(refC));
	    memcpy((char*)&altC,       buffer+10+5*j, sizeof(altC));
	    memcpy((char*)&cpgC,       buffer+12+5*j, sizeof(cpgC));

	    cout<<refC<<","<<altC<<":"<<(cpgC==1);
	    if(j < (sizePops -1 ))
	       cout<<"\t";

	}
	cout<<endl;


	
	if(0){
	cerr<<myFilezipped.good()<<endl;
	uint16_t chr;
	uint32_t coordinate ;
	char ref;
	char alt;

	myFilezipped.read((char*)&chr,        sizeof(chr));
	myFilezipped.read((char*)&coordinate, sizeof(coordinate));
	myFilezipped.read((char*)&ref,        sizeof(ref));
	myFilezipped.read((char*)&alt,        sizeof(alt));

	cout<<chrKnown[chr]<<"\t";
	cout<<coordinate<<"\t";
	cout<<"NACGT"[ref]<<","<<"NACGT"[alt]<<"\t";

	for(unsigned j=0;j<(sizePops);j++){
	    short refC;
	    short altC;
	    char  cpgC;

	    myFilezipped.read((char*)&refC,        sizeof(refC));
	    myFilezipped.read((char*)&altC,        sizeof(altC));
	    myFilezipped.read((char*)&cpgC,        sizeof(cpgC));

	    cout<<refC<<","<<altC<<":"<<(cpgC==1);
	    if(j < (sizePops -1 ))
	       cout<<"\t";

	}
	cout<<endl;
	cerr<<myFilezipped.good()<<endl;
	}
	//return 1;
	totalRecords++;

    }


    myFilezipped.close();
    cerr<<"Program "<<argv[0]<<" decompressed "<<totalRecords<<" terminated gracefully";

    return 0;
}

