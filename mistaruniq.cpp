/*
 * mistaruniq
 * Date: Apr-23-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

// #define DEBUG

#include "utils.h"
#include "MistarParser.h"
//#include "mistarunion.h"
#include "mistarOperations.h"

using namespace std;



int main (int argc, char *argv[]) {
    

    if(argc < 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr<<"usage: "<<argv[0]<<" [mistar file 1] [mistar file 2] ...\nwill print the unique union of the mistar files to stdout. Used to merge the same files from different filters"<<endl;
	return 1;       
    }

    




    vector<MistarParser * > vectorOfMP;
    for(int i=1;i<(argc);i++){ 
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
    cout<<"#MISTARUNIQ:"<<endl;
    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	cout<<"#UNIQFILE#"<<(i+1)<<endl;
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
	      coordCurrent,
	      true); //only print first

    bool atLeastOneHasData=true;///should be true

    bool stayLoop=true;
    
    //check if all have same pops
    vector<string> first=vector<string>(*vectorOfMP[0]->getPopulationsNames());
    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	vector<string> temp=vector<string>(*vectorOfMP[i]->getPopulationsNames());
	if(first !=  temp ){
	    cerr<<"Error file argument #"<<i<<" has different populations "<<vectorToString(first)<<" vs "<<vectorToString(temp)<<endl;
	    return 1;
	}
    }

    while(stayLoop){
	if(!atLeastOneHasData ){
	    stayLoop=false;
	    break;
	}

#ifdef DEBUG
	cerr<<"coordCurrent "<<coordCurrent<<endl;
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
	    string chrcheck = "";  
	    char refAllele  = '\0'; 
	    sanityCheck(vectorOfMP,
			hasData,
			hasCoordinate,
			vecAlleleRecords,
			chr1,
			coordCurrent,
			chrcheck,
			refAllele);
	    // printAllele(vectorOfMP,
	    // 		hasData,
	    // 		hasCoordinate,
	    // 		popSizePerFile,
	    // 		vecAlleleRecords,
	    // 		chr1,
	    // 		coordCurrent);
	    bool printedOne=false;
	    AlleleRecords arToCheck;
	    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
		if(hasData[i] && hasCoordinate[i]){
		   if(!printedOne){
		       cout<<(*vecAlleleRecords[i])<<endl;		       
		       printedOne=true;
		       arToCheck=(*vecAlleleRecords[i]);
		   }else{
		       if( !(arToCheck == (*vecAlleleRecords[i]) ) ){ //values should be equal 
			   cerr<<"Discrepency in data at coordinate "<<coordCurrent<<endl;
			   return 1;
		       }
		   }
		}
	    }

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
	     bool needToSetCoord=true;
	     for(unsigned int i=0;i<vectorOfMP.size();i++){ 

		if(hasData[i]  ){		    
		    if(needToSetCoord){
			coordCurrent  = vecAlleleRecords[i]->coordinate;
			needToSetCoord=false;
		    }else{
			coordCurrent  = min(coordCurrent,vecAlleleRecords[i]->coordinate);
		    }
		    // if(i==0){
		    // 	coordCurrent  = vecAlleleRecords[i]->coordinate;
		    // }else{
		    // 	coordCurrent  = min(coordCurrent,vecAlleleRecords[i]->coordinate);
		    // }	
		}
	     }

	    continue;

	}else{
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

