/*
 * replaceAncestor
 * Date: May-13-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>
#include <fstream>

// #define DEBUG

#include "utils.h"
#include "MistarParser.h"

using namespace std;

// inline unsigned int repositionCoordinate(const bool hasData1,const bool hasData2,const unsigned int record1Coordinate,const unsigned int record2Coordinate){

//     if( hasData1 && hasData2)
// 	return min( record1Coordinate,
// 		    record2Coordinate );
    
    
//     if(!hasData1 && hasData2)
// 	return record2Coordinate;
//     if( hasData1 && !hasData2)
// 	return record1Coordinate;

//     if( !hasData1 && !hasData2)
// 	return 0;
//     //this should never be reached
//     return 0;
// }


// //print an empty first record
// void printFirstEmpty(const unsigned int nonPop1, const unsigned int nonPop2, const AlleleRecords * record1,const AlleleRecords * record2){
//     cout<<record2->chr<<"\t"<<record2->coordinate<<"\t"<<record2->ref<<","<<record2->alt<<"\t";
//     //root
//     for(unsigned int i=0;i<1;i++){
// 	cout<<(*record2->vectorAlleles)[i];
// 	// stringify( ( (*record2->vectorAlleles)[i] ).refCount)<<","<<
// 	// stringify( ( (*record2->vectorAlleles)[i] ).altCount)<<":"<<
// 	// stringify( ( (*record2->vectorAlleles)[i] ).isCpg);
// 	// if( i!= (nonPop2-1) )
// 	cout<<"\t";
//     }


//     for(unsigned int i=1;i<nonPop1;i++){
// 	cout<<"0,0:0"<<"\t";
//     }

//     for(unsigned int i=1;i<nonPop2;i++){
// 	cout<<(*record2->vectorAlleles)[i];
// 	// stringify( ( (*record2->vectorAlleles)[i] ).refCount)<<","<<
// 	// stringify( ( (*record2->vectorAlleles)[i] ).altCount)<<":"<<
// 	// stringify( ( (*record2->vectorAlleles)[i] ).isCpg);
// 	if( i!= (nonPop2-1) )
// 	    cout<<"\t";
//     }
//     cout<<endl;
// }

// void printSecondEmpty(const unsigned int nonPop1, const unsigned int nonPop2, const AlleleRecords * record1,const AlleleRecords * record2){
//     cout<<record1->chr<<"\t"<<record1->coordinate<<"\t"<<record1->ref<<","<<record1->alt<<"\t";

//     for(unsigned int i=0;i<nonPop1;i++){
// 	cout<<(*record1->vectorAlleles)[i];
// 	// stringify( ( (*record1->vectorAlleles)[i] ).refCount)<<","<<
// 	// stringify( ( (*record1->vectorAlleles)[i] ).altCount)<<":"<<
// 	// stringify( ( (*record1->vectorAlleles)[i] ).isCpg);
// 	cout<<"\t";
//     }

//     for(unsigned int i=1;i<nonPop2;i++){
// 	cout<<"0,0:0";
// 	if( i!= (nonPop2-1) )
// 	    cout<<"\t";
//     }
//     cout<<endl;
// }


int main (int argc, char *argv[]) {
    
    // if(argc != 2 ){
    // cerr<<"usage: "<<argv[0]<<" [mistar file 1] [mistar file 2] \nwill print to stdout"<<endl;
    // return 1;
    // }

    if(argc != 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr<<"usage: "<<argv[0]<<" [mistar file 1] [mistar file 2]\nwill print the first mistar file but with the ancestral information from the second one to stdout"<<endl;
	return 1;
    }

    MistarParser mp1 (argv[1]);
    MistarParser mp2 (argv[2]);

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    cout<<"#MISTAR"<<endl;    	
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;
    cout<<"#REPLACEANCESTOR:"<<endl;
    //for(unsigned int i=0;i<vectorOfMP.size();i++){ 
    


    cout<<"#REPLACEANCESTOR#"<<(1)<<endl;
    cout<<""<<mp1.getHeader("#\t")<<"\n";
    cout<<"#REPLACEANCESTOR#"<<(2)<<endl;
    cout<<""<<mp2.getHeader("#\t")<<"\n";



    cout<<"#chr\tcoord\tREF,ALT\troot\tanc\t";

    //cout<<vectorToString( *(mp1.getPopulationsNames())," ")<<" "<<vectorToString( *(mp2.getPopulationsNames())," ")<<endl;
    bool hasData1 = mp1.hasData();
    bool hasData2 = mp2.hasData();

    unsigned int nonPop1 = mp1.getPopulationsNames()->size();
    //unsigned int nonPop2 = mp2.getPopulationsNames()->size();
    

    //printing first
    for(unsigned int i=2;i<nonPop1;i++){
	cout<<mp1.getPopulationsNames()->at(i);
	if( i!= (nonPop1-1) ){
	    cout<<"\t";
	}
    }
    cout<<endl;;
    // for(unsigned int i=1;i<nonPop2;i++){
    // 	cout<<mp2.getPopulationsNames()->at(i);
    // 	if(i!=(nonPop2-1)){
    // 	    cout<<"\t";
    // 	}
    // }


    if(!hasData1){
	cerr<<"File #1 does not have any data, exiting"<<endl;
	return 1;
    }
    if(!hasData2){
	cerr<<"File #2 does not have any data, exiting"<<endl;
	return 1;
    }

    AlleleRecords * record1 = mp1.getData();
    AlleleRecords * record2 = mp2.getData();

    if(record1->chr != record2->chr ){
	cerr<<"Chromosomes differ between "<<*record1<<" and "<<*record2<<endl;
	return 1;
    }

    // unsigned int coordCurrent=min( record1->coordinate,
    // 				   record2->coordinate );
    bool stayLoop=true;


    while(stayLoop){
	if(!hasData1  ){
	    stayLoop=false;
	    break;
	}


	//	return 1;

	//second one has data
	if(hasData2 &&
	   record1->coordinate == record2->coordinate ){
	    // return 1;
	    if(record1->chr != record2->chr ){
		cerr<<"Chromosomes differ between "<<(*record1)<<" and "<<(*record2)<<endl;
		return 1;
	    }


	    if(record1->ref != record2->ref ){
		cerr<<"The reference allele differs between "<<(*record1)<<" and "<<(*record2)<<endl;
		return 1;
	    }


	    
	    char newAlt='N';
	    // int rootLocated=0; //o=unknown,1=record1,2=record2,3=both should agree
	    // return 1;

	    if(record1->alt == record2->alt ){ //agree
		// rootLocated=3;
		newAlt=record1->alt;
		goto printnewrecord;
	    }
	    

	    if( record1->alt == 'N' && isResolvedDNA(record2->alt) ){
		// rootLocated=2;
		newAlt=record2->alt;
		goto printnewrecord;	
	    }

	    if(isResolvedDNA(record1->alt) && record2->alt == 'N' ){
		// rootLocated=1;
		newAlt=record1->alt;
		goto printnewrecord;	
	    }


	    //check for diff alt allele
	    if(record1->alt != record2->alt &&
	       isResolvedDNA(record1->alt)  &&
	       isResolvedDNA(record2->alt) ){ 

		//check if the alternative is due to the root or other populations
		int sumAlt=0;
		for(unsigned int i=0;i<2;i++){
		    sumAlt+=(*record2->vectorAlleles)[i].getAltCount();
		}

		
		if(sumAlt!=0) //the alt is either the root or anc
		    goto seekdata;//cannot reconcile alternative alleles

		//else safe ignore the second alternative since the second does not have the alt and pick the first as the alt
		newAlt=record1->alt;
		goto printnewrecord;//cannot reconcile alternative alleles
				    
		
	    }else{
		cerr<<"Wrong state between "<<(*record1)<<" and "<<(*record2)<<endl;
	    }

	    
	printnewrecord:
	    cout<<record1->chr<<"\t"<<record1->coordinate<<"\t"; //<<record2->coordinate<<"\t";
	    cout<<record1->ref<<",";
	    cout<<newAlt<<"\t";
	    
	    for(unsigned int i=0;i<2;i++){
		cout<<(*record2->vectorAlleles)[i];
		cout<<"\t";
	    }
	    
	    for(unsigned int i=2;i<record1->vectorAlleles->size();i++){
		cout<<(*record1->vectorAlleles)[i];
		if( i!= (record1->vectorAlleles->size()-1) )
		    cout<<"\t";
	    }
	    cout<<endl;

	    
	}
	else{//no record in the second one
	    
	    //check second record
 	    if(hasData2){
		//file 2 is behind, need to increase the 
		if( record1->coordinate > record2->coordinate ){
		    hasData2 = mp2.hasData();
		    if(hasData2){
			record2 = mp2.getData();
		    }	    
		    continue;//next iteration
		}
		    
 	    }
	    cout<<record1->chr<<"\t"<<record1->coordinate<<"\t"; //<<record2->coordinate<<"\t";
	    cout<<record1->ref<<",";
	    // cerr<<record1->chr<<"\t"<<record1->coordinate<<"\t"; //<<record2->coordinate<<"\t";
	    // cerr<<record1->ref<<",";
	    
	    // return 1;


	    int sumAlt=0;
	    for(unsigned int i=2;i<record1->vectorAlleles->size();i++){
		sumAlt+=(*record1->vectorAlleles)[i].getAltCount();
	    }
	    
	    if(sumAlt==0){ //the only alternative allele was due to the root
		cout<<"N\t";		
	    }else{
		cout<<record1->alt<<"\t";
	    }
	    SingleAllele s1;
	    cout<<s1<<"\t"<<s1<<"\t";

	    for(unsigned int i=2;i<record1->vectorAlleles->size();i++){
		cout<<(*record1->vectorAlleles)[i];
		if( i!= (record1->vectorAlleles->size()-1) )
		    cout<<"\t";
	    }
	    
	    cout<<endl;
	}

    seekdata:	
	if(hasData1){
	    hasData1 = mp1.hasData();
	    if(hasData1){
		record1 = mp1.getData();
	    }
	}



    }//end stayloop


    cerr<<"Program finished gracefully"<<endl;

    return 0;
}




	//repositioning the second one

// #ifdef DEBUG
// 	cerr<<"coordCurrent "<<coordCurrent<<endl;
// #endif

// 	//overlap between data
// 	if(coordCurrent == record1->coordinate &&
// 	   coordCurrent == record2->coordinate ){


// #ifdef DEBUG
// 	    cerr<<"enter "<<coordCurrent<<endl;
// #endif


// 	    if(record1->chr != record2->chr ){
// 		cerr<<"Chromosomes differ between "<<alleleRecordsAsString(*record1)<<" and "<<alleleRecordsAsString(*record2)<<endl;
// 		return 1;
// 	    }


// 	    if(record1->ref != record2->ref ){
// 		cerr<<"The reference allele differs between "<<alleleRecordsAsString(*record1)<<" and "<<alleleRecordsAsString(*record2)<<endl;
// 		return 1;
// 	    }




// 	//     char newAlt;
// 	//     int rootLocated=0; //o=unknown,1=record1,2=record2,3=both should agree

// 	//     if(record1->alt == record2->alt ){ //agree
// 	// 	rootLocated=3;
// 	// 	newAlt=record1->alt;
// 	// 	goto printnewrecord;
// 	//     }

// 	//     if( record1->alt == 'N' || isResolvedDNA(record2->alt) ){
// 	// 	rootLocated=2;
// 	// 	newAlt=record2->alt;
// 	// 	goto printnewrecord;	
// 	//     }

// 	//     if(isResolvedDNA(record1->alt) || record2->alt == 'N' ){
// 	// 	rootLocated=1;
// 	// 	newAlt=record1->alt;
// 	// 	goto printnewrecord;	
// 	//     }

// 	//     //cannot reconcile alternative alleles
// 	//     goto seekdata;


// 	// printnewrecord:
// 	//     cout<<record1->chr<<"\t"<<record1->coordinate<<"\t"<<record1->ref<<","<<newAlt<<"\t";

// 	//     //print root
// 	//     switch(rootLocated){
// 	//     case 1:
// 	// 	for(unsigned int i=0;i<1;i++){
// 	// 	    cout<<(*record1->vectorAlleles)[i];
// 	// 	    // stringify( ( (*record1->vectorAlleles)[i] ).refCount)<<","<<
// 	// 	    // stringify( ( (*record1->vectorAlleles)[i] ).altCount)<<":"<<
// 	// 	    // stringify(( (*record1->vectorAlleles)[i] ).isCpg);
// 	// 	    cout<<"\t";
// 	// 	}

// 	// 	break;

// 	//     case 2:
// 	// 	for(unsigned int i=0;i<1;i++){
// 	// 	    cout<< (*record2->vectorAlleles)[i];
// 	// 	    // stringify( ( (*record2->vectorAlleles)[i] ).refCount)<<","<<
// 	// 	    // stringify( ( (*record2->vectorAlleles)[i] ).altCount)<<":"<<
// 	// 	    // stringify(( (*record2->vectorAlleles)[i] ).isCpg);
// 	// 	    cout<<"\t";
// 	// 	}

// 	// 	break;

// 	//     case 3:

// 	// 	if((*record1->vectorAlleles)[0] != (*record2->vectorAlleles)[0]){
// 	// 	    cerr<<"Both records disagree on the ancestral allele information at "<<record1->coordinate<<endl;
// 	// 	    return 1;
// 	// 	}



// 	// 	//agree on the root
// 	// 	for(unsigned int i=0;i<1;i++){
// 	// 	    cout<<(*record1->vectorAlleles)[i];
// 	// 	    // stringify( ( (*record1->vectorAlleles)[i] ).refCount)<<","<<
// 	// 	    // stringify( ( (*record1->vectorAlleles)[i] ).altCount)<<":"<<
// 	// 	    // stringify(( (*record1->vectorAlleles)[i] ).isCpg);
// 	// 	    cout<<"\t";
// 	// 	}


// 	// 	break;

// 	//     default:
// 	// 	cerr<<"Both records disagree on the ancestral allele at "<<record1->coordinate<<endl;
// 	// 	return 1;
// 	// 	break;	


// 	//     }

// 	//     for(unsigned int i=1;i<record1->vectorAlleles->size();i++){
// 	// 	cout<<(*record1->vectorAlleles)[i];
// 	// 	cout<<"\t";
// 	//     }

// 	//     for(unsigned int i=1;i<record2->vectorAlleles->size();i++){
// 	// 	cout<<(*record2->vectorAlleles)[i];
// 	// 	if( i!= (record2->vectorAlleles->size()-1) )
// 	// 	    cout<<"\t";
// 	//     }
// 	//     cout<<endl;

// 	seekdata:

// 	    if(hasData1){
// 		hasData1 = mp1.hasData();
// 		if(hasData1){
// 		    record1 = mp1.getData();
// 		}
// 	    }

// 	    if(hasData2){
// 		hasData2 = mp2.hasData();
// 		if(hasData2){
// 		    record2 = mp2.getData();
// 		}
// 	    }

// 	    coordCurrent=repositionCoordinate(hasData1,hasData2,record1->coordinate,record2->coordinate);
// 	    if(coordCurrent == 0){
// 		cerr<<"Program finished gracefully";
// 		return 0;
// 	    }





// 	}else{ //different coordinates

// #ifdef DEBUG
// 	    cerr<<"diff "<<coordCurrent<<" "<<hasData1<<" "<<hasData2<<" record 1 "<<record1->coordinate<< " record 2 "<<record2->coordinate<<endl;
// #endif







// 	    //if the counter is ahead of both of them and both still have data
// 	    if( ((coordCurrent > record1->coordinate) && hasData1) ||
// 		((coordCurrent > record2->coordinate) && hasData2)){ //running behind for one of them
// 		cerr<<"Internal error: Are the mistar files sorted by coordinate? at "<<coordCurrent<< " record 1 "<<record1->coordinate<< " record 2 "<<record2->coordinate<<endl;
// 		return 1;
// 	    }



// 	    //if they both have data
// 	    if(hasData1 && hasData2){
// 		//print record1, no data for 2
// 		if(coordCurrent == record1->coordinate &&
// 		   coordCurrent < record2->coordinate){
// #ifdef DEBUG
// 		    cerr<<"case1 "<<coordCurrent<<endl;
// #endif



// 		    printSecondEmpty(nonPop1,nonPop2,record1,record2);


// 		    if(hasData1){
// 			hasData1 = mp1.hasData();
// 			if(hasData1){
// 			    record1 = mp1.getData();
// 			}
// 		    }

// 		    coordCurrent=repositionCoordinate(hasData1,hasData2,record1->coordinate,record2->coordinate);
// 		    if(coordCurrent == 0){
// 			cerr<<"Program finished gracefully";
// 			return 0;
// 		    }


// 		    continue;

// 		}





// 		if(coordCurrent < record1->coordinate &&
// 		   coordCurrent == record2->coordinate){
// #ifdef DEBUG
// 		    cerr<<"case2 "<<coordCurrent<<endl;
// #endif


// 		    printFirstEmpty(nonPop1,nonPop2,record1,record2);

// 		    if(hasData2){
// 			hasData2 = mp2.hasData();
// 			if(hasData2){
// 			    record2 = mp2.getData();
// 			}
// 		    }

// 		    coordCurrent=repositionCoordinate(hasData1,hasData2,record1->coordinate,record2->coordinate);
// 		    if(coordCurrent == 0){
// 			cerr<<"Program finished gracefully";
// 			return 0;
// 		    }

// 		    continue;
// 		}
// 	    }else{

// 		//one of them does not have data
// 		if(!hasData1 && hasData2){
// 		    printFirstEmpty(nonPop1,nonPop2,record1,record2);
// 		    hasData2 = mp2.hasData();
// 		    if(hasData2){
// 			record2 = mp2.getData();
// 		    }

// 		}

// 		if( hasData1 && !hasData2){	
// 		    printSecondEmpty(nonPop1,nonPop2,record1,record2);

// 		    hasData1 = mp1.hasData();
// 		    if(hasData1){
// 			record1 = mp1.getData();
// 		    }

// 		}


// 		coordCurrent=repositionCoordinate(hasData1,hasData2,record1->coordinate,record2->coordinate);


// 	    }



// 	}//end not correct coord

