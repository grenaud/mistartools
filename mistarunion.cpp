/*
 * mistarunion
 * Date: Jan-25-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

//#define DEBUG

#include "utils.h"
#include "MistarParser.h"
//#include "mistarunion.h"
#include "mistarOperations.h"

using namespace std;



int main (int argc, char *argv[]) {
    bool force=false;    

    if(argc < 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr<<"usage: "<<argv[0]<<" <options> [mistar file 1] [mistar file 2] ...\nwill print the union of the mistar files to stdout, it will skip triallelic sites"<<endl<<endl<<
	    "\tOptions:\n"<<
	    "\t\t-f\t\tForce ignore of error like different reference alleles"<<
	    endl;
	return 1;       
    }

    
    



    vector<MistarParser * > vectorOfMP;
    for(int i=1;i<(argc);i++){ 

	if(i==1 && string(argv[i]) == "-f"){
	    force=true;
	    continue;
	}

	MistarParser * mp = new MistarParser(argv[i]);
	vectorOfMP.push_back(mp);
    }    


    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    cout<<"#MISTAR"<<endl;    	
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;
    cout<<"#MISTARUNION:"<<endl;
    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	cout<<"#UNIONFILE#"<<(i+1)<<endl;
	cout<<""<<vectorOfMP[i]->getHeader("#\t")<<"\n";
    }
    

    vector<bool> hasData;
    vector<int> popSizePerFile;
    vector<AlleleRecords *> vecAlleleRecords;
    string chr1;
    unsigned int coordCurrent;
    
    initFiles(vectorOfMP,
	      //atLeastOneHasData,
	      hasData,
	      popSizePerFile,
	      vecAlleleRecords,
	      chr1,
	      coordCurrent);

    bool atLeastOneHasData=true;///should be true

    bool stayLoop=true;


    while(stayLoop){
	if(!atLeastOneHasData ){
	    stayLoop=false;
	    break;
	}

#ifdef DEBUG
	cerr<<"coordCurrent "<<chr1<<"\t"<<coordCurrent<<endl;
#endif

	vector<bool> hasCoordinate (vectorOfMP.size(),false);
	bool atLeastOneHasCoordinate=false;
	for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	    if(hasData[i]){
		if(coordCurrent  == vecAlleleRecords[i]->coordinate){
		    hasCoordinate[i]=true;
		    atLeastOneHasCoordinate=true;
		}
	    }
	}
	// cout<<vectorToString(hasData,"-")<<endl;

	
	//we print
	if(atLeastOneHasCoordinate){
	    
	    printAllele(vectorOfMP,
			hasData,
			hasCoordinate,
			popSizePerFile,
			vecAlleleRecords,
			chr1,
			coordCurrent,	    
			force);


	    // 	seekdata:
	     atLeastOneHasData=false;

	     for(unsigned int i=0;i<vectorOfMP.size();i++){ 
		 
		 if(hasData[i] ){


		     //only get data from those with the coordinate
		     if(hasCoordinate[i]){
			 hasData[i]  =  vectorOfMP[i]->hasData();
			 if(hasData[i]){
			     atLeastOneHasData=true;
			     vecAlleleRecords[i] = vectorOfMP[i]->getData() ;
			 }
		     }else{
			 atLeastOneHasData=true;//still one with data
		     }

		}
	     }


	     // cout<<"seek "<<vectorToString(hasData,"-")<<"\t"<<atLeastOneHasData<<endl;

	     if(!atLeastOneHasData)
		 continue;

	     //coordCurrent=0;
	     //Try to find the record that is the most behind
	     bool needToSetCoord=true;
	     for(unsigned int i=0;i<vectorOfMP.size();i++){ 

		if(hasData[i]  ){		    
		    if(needToSetCoord){

			chr1          = vecAlleleRecords[i]->chr;
			coordCurrent  = vecAlleleRecords[i]->coordinate;

			needToSetCoord=false;
		    }else{
			if(chr1 == vecAlleleRecords[i]->chr){
			    coordCurrent  = min(coordCurrent,vecAlleleRecords[i]->coordinate);
			}else{

			    int chrcmp = compare2Chrs(chr1, vecAlleleRecords[i]->chr);
			    if(chrcmp == -1){// the current record is ahead, do nothing
				
				
			    }else{
				if(chrcmp== 1){ //the current record is behind chromosome-wise, reposition there.
				    chr1          = vecAlleleRecords[i]->chr;
				    coordCurrent  = vecAlleleRecords[i]->coordinate;

				}else{
				    cerr<<"Invalid state"<<endl;
				    return 1;  				    
				}
			    }
			}
		    }
		    // if(i==0){
		    // 	coordCurrent  = vecAlleleRecords[i]->coordinate;
		    // }else{
		    // 	coordCurrent  = min(coordCurrent,vecAlleleRecords[i]->coordinate);
		    // }	
		}
	     }

	    continue;

	}else{//not at least one has coordinate
	     cerr<<"Invalid state"<<endl;
	     return 1;  
	}

	
    }


    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	//MistarParser * mp = new MistarParser(argv[i]);
	delete(vectorOfMP[i]);
    }    


    cerr<<"Program finished gracefully"<<endl;



    return 0;
}

