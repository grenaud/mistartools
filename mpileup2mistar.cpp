/*
 * mpileup2mistar
 * Date: Jun-21-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>

#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {

    if(argc != 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
        cerr << "This program creates a mistar matrix using the output for a single sample:\nsamtools mpileup  -f ref.fa filein.bam\n\n\nUsage: "<<argv[0]<<" [samtools mpileup] [sample name]"<<endl;
	return 1;        
    }

    vector<unsigned int> char2index (99,90);
    for(unsigned int i=0;i<90;i++){
	if(i==65)
	    char2index[i] = 0;
	if(i==67)
	    char2index[i] = 1;
	if(i==71)
	    char2index[i] = 2;
	if(i==84)
	    char2index[i] = 3;
    }
    


    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

   string filename   = string(argv[argc-2]);
   string samplename = string(argv[argc-1]);


    cout<<"#MISTAR"<<endl;    	
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;
	
    cout<<"#chr"<<"\t"<<"coord"<<"\t"<<"REF,ALT\troot\tanc\t"<<samplename<<endl;


   string line;
   igzstream myFile;

   myFile.open(filename.c_str(), ios::in);

   if (myFile.good()){
       while ( getline (myFile,line)){
	   vector<string> tempfield   = allTokens( line,'\t');
	   vector<unsigned int> counterAllele (0,4);

	   //0          1       2       3       4               5
	   //chr10	6	A	9	......,,.	D:88CB16=
	   char ref               = tempfield[2][0];
	   unsigned int refIndex  = char2index[ ref ];
	   if(!isResolvedDNA(ref))//only where the reference is A,C,G,T
	       continue;

	   //       chr                   coord
	   //cout<<tempfield[0]<<"\t"<< tempfield[1] <<"\t";
	   
	   for(unsigned i=0;i<tempfield[4].size();i++){
	       switch(tempfield[4][i]){
	       case '.':  //increase ref
		   counterAllele[refIndex]++;
		   break;

	       case ',': //increase ref
		   counterAllele[refIndex]++;
		   break;
		   
		   //beginning of read
	       case '^': // ignore ! in next : ^! 
		   i++; //skipping the !
		   break;

		   //end read
	       case '$': //ignore 
		   break;

		   
	       case '-': //insertions
		   break;

	       case '+': //insertions
		   break;
		   
	       default:
		   counterAllele[ char2index[ tempfield[4][i] ]  ]++;

		   break;
	       }
	   }
	   
       }
     myFile.close();
   }else{
       cerr << "Unable to open file "<<filename<<endl;
       return 1;
    }
    return 0;
}

