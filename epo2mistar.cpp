/*
 * epo2mistar
 * Date: Mar-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>

#include "utils.h"
#include "SingleAllele.h"

using namespace std;


inline int bp2Index(const char toCheck){
    if(toCheck == 'A')
	return 1;
    if(toCheck == 'C')
	return 2;
    if(toCheck == 'G')
	return 3;
    if(toCheck == 'T')
	return 4;

    return 0;
}


//Creates SingleAllele object from the information 
SingleAllele saFromInformation (char alleleEPO,
				char refAllele,
				char altAllele,
				bool ancestralCpG
				,string & line
				){

    if(refAllele == 'N' && altAllele == 'N' ){
	SingleAllele toReturn (0,
			       0,
			       ancestralCpG);
	return toReturn;
    }

    if(refAllele != 'N' && altAllele == 'N' ){
	// if(alleleEPO != refAllele){
	//     cerr<<"ERROR 3"<<line<<endl;
	//     exit(1);
	// }
	SingleAllele toReturn (alleleEPO == refAllele,
			       0,
			       ancestralCpG);
	return toReturn;
    }

    if(refAllele == 'N' && altAllele != 'N' ){
	cerr<<"ERROR 4 with line "<<line<<endl;
	exit(1);
    }

    if(refAllele != 'N' && altAllele != 'N' ){
	SingleAllele toReturn (alleleEPO == refAllele,
			       alleleEPO == altAllele,
			       ancestralCpG);
	return toReturn;
    }

    
    cerr<<"ERROR 5 with line "<<line<<endl;
    exit(1);
    
    SingleAllele toReturn (0,
			   0,
			   ancestralCpG);
    return toReturn;
}

int main (int argc, char *argv[]) {
    if(argc != 2 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr<<"usage: "<<argv[0]<<" [parsed epo file]\nwill print the mistar files to stdout"<<endl;
	return 1;       
    }

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    cout<<"#MISTAR"<<endl;    	
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;
	
    cout<<"#chr"<<"\t"<<"coord"<<"\t"<<"REF,ALT\troot\tanc\thref\tchimpAnc\tchimp\tGorAnc\tGor\tOrangAnc\tOrang"<<endl;

    string line;
    igzstream myFile;
    string filename = string(argv[1]);
    myFile.open(filename.c_str(), ios::in);
    
    unsigned int totalRec=0;
    unsigned int writtenRec=0;

    
    if (myFile.good()){
	while ( getline (myFile,line)){
	    totalRec++;

	    vector<string> fields=allTokens(line,'\t');
	    bool   ancestralCpG=false;

	    //filtering out triallelic sites
	    bool  altAllelesFlags [5] = {false,false,false,false,false};  //will be true once we encounter certain alt alleles
	    for(unsigned i=2;i<=8;i++){
		altAllelesFlags[ bp2Index( fields[i][0] ) ]=true;
	    }

	    int numberOfAltAlleles=0;
	    for(int i=1;i<=4;i++){
		if(altAllelesFlags[i])
		    numberOfAltAlleles++;
	    }

	    if(numberOfAltAlleles>2)
		continue;



	    if(fields[9][0] == '1') //Ancestral is CpG
		ancestralCpG=true;
	    
	    char refAllele =  fields[2][0]  ;

	    if(refAllele == 'N'){//we should not trust the alignment if the refence is not defined
		continue;
	    }

	    int  refAlleleIdx= bp2Index(refAllele);
	    char altAllele = 'N';

	    for(int i=1;i<=4;i++){
		if(altAllelesFlags[i] && i!=refAlleleIdx)
		    altAllele="NACGT"[i];
	    }
	    int  altAlleleIdx= bp2Index(altAlleleIdx);

	    

	    cout<<fields[0]<<"\t"<<fields[1]<<"\t"<<refAllele<<","<<altAllele<<"\t";

	    vector<SingleAllele> toprintpop;

	    //root (chimp)
	    toprintpop.push_back( saFromInformation(fields[4][0],refAllele,altAllele,ancestralCpG,line ) );

	    //anc
	    toprintpop.push_back( saFromInformation(fields[3][0],refAllele,altAllele,ancestralCpG,line ) );
	    
	    //href
	    toprintpop.push_back( saFromInformation(fields[2][0],refAllele,altAllele,ancestralCpG,line ) );

	    //anc
	    toprintpop.push_back( saFromInformation(fields[3][0],refAllele,altAllele,ancestralCpG,line ) );
	    
	    //root (chimp)
	    toprintpop.push_back( saFromInformation(fields[4][0],refAllele,altAllele,ancestralCpG,line ) );

	    //remaining
	    for(unsigned i=5;i<=8;i++){
		toprintpop.push_back( saFromInformation(fields[i][0],refAllele,altAllele,ancestralCpG,line) );
	    }
	    
	    writtenRec++;

	    cout<<vectorToString(toprintpop,"\t")<<endl;
	    

	}
	myFile.close();
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	return 1;
    }

    cerr<<"Program "<<argv[0]<<" looked at  "<<totalRec<<" records, wrote "<<writtenRec<<" terminated gracefully"<<endl;

    return 0;
}

