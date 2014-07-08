/*
 * parser
 * Date: Apr-23-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>

#include "utils.h"
#include "MSParser.h"
#include "MSobject.h"
#include "NjTree.h"


using namespace std;

int main (int argc, char *argv[]) {

    if(argc == 1){
        cerr << "This program creates a simple neighbor joining tree using ms output\n\n!!! Please note that the first individual will be used as root !!!\n\nUsage: "<<argv[0]<<" [ms file#1] [ms file#2] ... "<<endl;
       return 1;        
    }



    
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    int numberOfIndividuals=-1;
    vector< vector<double> > distanceMatrix ; //= vector< vector<double> > (numberOfPopulations, vector<double>(numberOfPopulations));
    vector<string> namesToUse ; // = vector<string>(*names);


    for(int arg=1;arg<argc;arg++){
	cerr<<"opening "<<string(argv[arg]) <<endl;

	MSParser msp (    string(argv[arg]) );
	for(int numRec=0;numRec<msp.numberOfRecords();numRec++){

	    //cerr<<"Looking at record# "<<numRec<<endl;

	    // if(msp.numberOfRecords() != 1){
	    //     cerr<<"error wrong # records "<<string(argv[arg]) <<endl;
	    //     return 1;
	    // }
	    const MSobject * mso = msp.getMSObj(numRec);
	    if(numberOfIndividuals==-1){
		numberOfIndividuals=mso->getNumberIndividuals();
		distanceMatrix = vector< vector<double> > (numberOfIndividuals, vector<double>(numberOfIndividuals,0.0));
		namesToUse.push_back("root");
		for(int i=1;i<numberOfIndividuals;i++){
		    namesToUse.push_back("i"+stringify(i));
		}
	    }else{
		if(mso->getNumberIndividuals() != numberOfIndividuals){
		    cerr<<"error wrong # individuals "<<string(argv[arg]) <<" found "<<mso->getNumberIndividuals()<<" and not "<<numberOfIndividuals<<endl;
		    return 1;
		}	    
	    }

	    // cout<<".";
	
	    vector<const bool *> indAllele;
	    for(int i=1;i<=mso->getNumberIndividuals();i++){
		indAllele.push_back(mso->getBoolArrayIndividual(i));
	    }
	    // cout<<".";
	
	    for(unsigned int i=0;i<mso->getNumberSegSites();i++){

		for(int j=0;j<numberOfIndividuals;j++){
		    for(int k=0;k<numberOfIndividuals;k++){
			if(j!=k){
			    if(indAllele[j][i] != indAllele[k][i]){
				distanceMatrix[k][j]++;
				distanceMatrix[j][k]++;
			    }
			}
		    }
		}
	    }
	    // cout<<"..done"<<endl;
	}
	    
    }



    //cout<<"\t"<<vectorToString(namesToUse,"\t")<<endl;
    if(0)
    for(int j=0;j<numberOfIndividuals;j++){
     	for(int k=0;k<numberOfIndividuals;k++){	   
	    // cout<<namesToUse[j]<<"\t"<<namesToUse[k]<<"\t"<<distanceMatrix[j][k]<<endl;
	}
     	//cout<<endl;
	//cout<<namesToUse[j]<<"\t"<<vectorToString(distanceMatrix[j],"\t")<<endl;
    }

    Tree * mytree=neighborJoinFromDist(distanceMatrix,namesToUse,numberOfIndividuals);
    cout<<*mytree<<endl;


    return 0;
}

