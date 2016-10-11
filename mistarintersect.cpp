/*
 * mistarintersect
 * Date: Jan-25-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

//#define DEBUG

#include "utils.h"
#include "MistarParser.h"
#include "mistarOperations.h"

using namespace std;



int main (int argc, char *argv[]) {
    
    bool force=false;    

    if(argc < 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr<<"usage: "<<argv[0]<<" [mistar file 1] [mistar file 2] ...\nwill print the intersection of the mistar files to stdout, it will skip triallelic sites"<<endl;
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
    cout<<"#MISTARINTERSECT:"<<endl;
    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	cout<<"#INTERSECTFILE#"<<(i+1)<<endl;
	cout<<""<<vectorOfMP[i]->getHeader("#\t")<<"\n";
    }
    
    //    bool atLeastOneHasData;
    vector<bool> hasData;
    vector<int> popSizePerFile;
    vector<AlleleRecords *> vecAlleleRecords;
    string chr1;
    unsigned int coordCurrent;

    
    initFiles(vectorOfMP,
	      // atLeastOneHasData,
	      hasData,
	      popSizePerFile,
	      vecAlleleRecords,
	      chr1,
	      coordCurrent);

    vector<bool>  hasCoordinate (vectorOfMP.size(),true);//dummy value


    bool stayLoop=true;


    while(stayLoop){

#ifdef DEBUG
	cerr<<"coordCurrent "<<chr1<<":"<<coordCurrent<<endl;
#endif


	bool allHaveCoordinate=true;

	for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	    if(hasData[i]){		
		if( (chr1          != vecAlleleRecords[i]->chr       ) ||
		    (coordCurrent  != vecAlleleRecords[i]->coordinate) ){
		    allHaveCoordinate=false;
		}
	    }else{
		stayLoop=false;
		break;
	    }
	}

	
	//we print
	if(allHaveCoordinate){
#ifdef DEBUG
	cerr<<"same coordCurrent "<<chr1<<":"<<coordCurrent<<endl;
#endif

	    
	    printAllele(vectorOfMP,
			hasData,
			hasCoordinate,
			popSizePerFile,
			vecAlleleRecords,
			chr1,
			coordCurrent,
			force);


	    // 	seekdata:
	    allHaveCoordinate=false;
	    
	    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
		 
		if(!hasData[i] ){
		    cerr<<"Invalid state"<<endl;
		    return 1;
		}
		hasData[i]  =  vectorOfMP[i]->hasData();
		if(hasData[i]){
		    vecAlleleRecords[i] = vectorOfMP[i]->getData() ;
		}else{
		    stayLoop=false;
		    break;
		}
		
	    }


	    //all have had getData called, we need to reposition to the maximum coord
	    bool needToSetCoord=true;
	    for(unsigned int i=0;i<vectorOfMP.size();i++){ 

#ifdef DEBUG
		cerr<<needToSetCoord<<"\t"<<i<<"\t"<<chr1<<":"<<coordCurrent<<endl;
#endif

		 if(needToSetCoord){
		     chr1          = vecAlleleRecords[i]->chr;
		     coordCurrent  = vecAlleleRecords[i]->coordinate;
		     needToSetCoord=false;
		 }else{
		     int chrcmp = compare2Chrs(chr1, vecAlleleRecords[i]->chr);

		     if(chrcmp != 0){//chromosomes are not equal

			 if(chrcmp == -1){ //chr1 < vecAlleleRecords[i]->chr
			     chr1          = vecAlleleRecords[i]->chr;
			     coordCurrent  = vecAlleleRecords[i]->coordinate;			     
			 }

			 if(chrcmp ==  1){ // chr1 > vecAlleleRecords[i]->chr
			      // chr1          = chr1;
			      // coordCurrent  = vecAlleleRecords[i]->coordinate;
			 }

			 

		     }else{
			 coordCurrent  = max(coordCurrent,vecAlleleRecords[i]->coordinate); //chromosomes are equal, jump to the max
		     }


		     // if( chr1          != vecAlleleRecords[i]->chr ){//need to skip to a new chr
		     // 	 if(coordCurrent  > vecAlleleRecords[i]->coordinate){//the diff chr is probably a new chr
		     // 	     chr1          = chr1;
		     // 	     coordCurrent  = vecAlleleRecords[i]->coordinate;
		     // 	 }
		     // }else{
		     // 	 if( chr1          == vecAlleleRecords[i]->chr ){
		     // 	     coordCurrent  = max(coordCurrent,vecAlleleRecords[i]->coordinate);
		     // 	 }else{//chr1 is greater than vecAlleleRecords[i]
		     // 	     //	we will reposition vecAlleleRecords[i] there 
		     // 	 }
			 
		     // }
		 }

	     }

	    continue;

	}else{
	     // cerr<<"Invalid state"<<endl;
	     // return 1;  
	    

#ifdef DEBUG
	    cerr<<"current "<<chr1<<"\t"<<coordCurrent<<endl;
#endif
	    
	    for(unsigned int i=0;i<vectorOfMP.size();i++){ 

#ifdef DEBUG
		cerr<<"coord["<<i<<"] "<<vecAlleleRecords[i]->chr<<"\t"<<vecAlleleRecords[i]->coordinate<<endl;
#endif

		//record [i] is ahead, re-position there
		if( ((chr1          == vecAlleleRecords[i]->chr) && 
		     (coordCurrent  < vecAlleleRecords[i]->coordinate)) ){ //overshot [i], repositioning there
		    chr1          = vecAlleleRecords[i]->chr       ;
		    coordCurrent  = vecAlleleRecords[i]->coordinate;
		    continue;
		}

		int chrcmp = compare2Chrs(chr1, vecAlleleRecords[i]->chr);

		//record [i] is ahead, re-position there
		if( chrcmp  ==  -1    ){// chr1 < vecAlleleRecords[i]->chr
		//if( chr1          != vecAlleleRecords[i]->chr         ){// chr1 < vecAlleleRecords[i]->chr

		    //if(coordCurrent  > vecAlleleRecords[i]->coordinate)  { //the different chr is probably a new chr
		    chr1          = vecAlleleRecords[i]->chr      ;
		    coordCurrent  = vecAlleleRecords[i]->coordinate;
		    continue;
		    //}

		}

		//if( (chr1         == vecAlleleRecords[i]->chr) &&
		//in sync
		if( (chrcmp        == 0) &&
		    (coordCurrent == vecAlleleRecords[i]->coordinate)){ //fine
		}
		    
		//if( ((chr1          != vecAlleleRecords[i]->chr) && (coordCurrent <  vecAlleleRecords[i]->coordinate)) 
		//record [i] is running behind 
		if( (chrcmp == 1)  // chr1 > vecAlleleRecords[i]->chr
		    ||
		    ((chr1          == vecAlleleRecords[i]->chr) && (coordCurrent >  vecAlleleRecords[i]->coordinate)) ){ 
		    hasData[i]  =  vectorOfMP[i]->hasData();
		    if(hasData[i]){
			vecAlleleRecords[i] = vectorOfMP[i]->getData() ;
		    }else{
			stayLoop=false;
			break;
		    }
		}
		
	    }

	    
	}//end different coord

	
    }//end main loop

    // finish:
    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	//MistarParser * mp = new MistarParser(argv[i]);
	delete(vectorOfMP[i]);
    }    


    cerr<<"Program finished gracefully"<<endl;



    return 0;
}

