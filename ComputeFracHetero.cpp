/*
 * ComputeDivergence
 * Date: Aug-20-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "ComputeFracHetero.h"

// #define DEBUG
// #define DEBUG2
// #define DEBUG3

// int computeFracHetero(string refereVCF,string refereVCFidx,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,GenomicRange grg,int minGQcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minMQcutoff,double minMapabilitycutoff,int bpForIndels){
//     return computeFracHetero( refereVCF, refereVCFidx, epoFile, epoFileidx, minimumNumberOfGoodSites, grg.getChrName(), grg.getStartCoord(), grg.getEndCoord(),  minGQcutoff, minCovcutoffREF, maxCovcutoffREF, minMQcutoff, minMapabilitycutoff, bpForIndels);

// }


int computeFracHetero(string refereVCF,string refereVCFidx,int minimumNumberOfGoodSites,GenomicRange grg,SetVCFFilters * filtersVCF, int bpForIndels){
		      //int minGQcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minMQcutoff,double minMapabilitycutoff,int bpForIndels){
    return computeFracHetero( refereVCF, refereVCFidx, minimumNumberOfGoodSites, grg.getChrName(), grg.getStartCoord(), grg.getEndCoord(), filtersVCF, bpForIndels);
			      //minGQcutoff, minCovcutoffREF, maxCovcutoffREF, minMQcutoff, minMapabilitycutoff, bpForIndels);

}


//int computeFracHetero( string refereVCF,string refereVCFidx,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,string chrName, int startChrCoord,  int endChrCoord,int minGQcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minMQcutoff,double minMapabilitycutoff,int bpForIndels){
int computeFracHetero(string refereVCF,string refereVCFidx,int minimumNumberOfGoodSites,string chrName, int startChrCoord,  int endChrCoord,SetVCFFilters * filtersVCF, int bpForIndels){
		       //int minGQcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minMQcutoff,double minMapabilitycutoff,int bpForIndels){

    int numberOfGoodSites=0;
	  

    //READING VCF DATA
    VCFreader vcfREF ( refereVCF.c_str(), refereVCFidx.c_str(), chrName, startChrCoord,endChrCoord,bpForIndels );
	    
    SimpleVCF * smvcfREF;

    bool lineLeftREF=vcfREF.hasData();
    // bool lineLeftEPO;

    if(lineLeftREF){ smvcfREF=vcfREF.getData();  } //calling the copy constructor...


    // //READING EPO DATA CHECK IF FLAG FAIL

    // ReadTabix rtEPO ( epoFile.c_str()  , epoFileidx.c_str()  , chrName, startChrCoord, endChrCoord );
    // string lineFromEPO;
    // lineLeftEPO=(rtEPO.readLine( lineFromEPO ));



	 

	    

    if(!( lineLeftREF  )){ 
	//continue; //need all lines
	return -1;
    }

    //    unsigned int coordCurrent=chrCoord;
    unsigned int coordCurrent=startChrCoord;

    
    int countHetero = 0;
    int countHomo   = 0;

    //DEBUG 
    // int rejectEPOValidREF=0;
    // int rejectEPOValidANC=0;

    // int rejectERROR_REF=0;
    // int rejectREF_genotypeUnknown=0;

	    

    bool stayLoop=true;
    while(stayLoop){
	//reached the end
	if( !lineLeftREF  ){
	    stayLoop=false;
	    break;
	}
			   
	// vector<string> fieldsEPO=allTokens(lineFromEPO,"\t");

	//current position in the VCF files
	unsigned int positionREF=smvcfREF->getPosition();
	// unsigned int positionEPO=string2uint(fieldsEPO[1]);
		    
#ifdef DEBUG2
	cout<<endl<<"######"<<endl;
	cout<<coordCurrent<<endl;
	cout<<positionREF<<endl;
	//	cout<<positionEPO<<endl;
#endif		


	//all at the same genomic coord
	if(coordCurrent == positionREF  ){ 
#ifdef DEBUG2
	    cout<<"entered"<<endl;
#endif
	    //increase for next iteration
	    coordCurrent++;






	    //////////////////////////////////////
	    //                                  //
	    //   BEGIN  VCF ACTUAL PROCESSING   //
	    //                                  //
	    //////////////////////////////////////

			
	    // //DETERMINING CPG
	    // bool isCpG=false;
	    // if(fieldsEPO[9]     == "0"){
	    // 	//TODO: CHECK previous alleles 
	    // 	//immediately to see if they were a C
			    
	    // }else{
	    // 	if(fieldsEPO[9] == "1"){
	    // 	    isCpG=true;
	    // 	}else{
	    // 	    cerr<<"ERROR: wrong CpG flag in EPO"<<endl;
	    // 	    return 1;
	    // 	}
	    // }

		

	    /////////////////////////////////////////////////////////////////
	    //           FILTERING INDELS AND UNDEFINED REF ALLELE         //
	    /////////////////////////////////////////////////////////////////

	    //passedFilters( SimpleVCF smvcf,int minCovcutoff,int maxCovcutoff,double minMapabilitycutoff,int minMQcutoff,int minGQcutoff);
	    //bool passingFilters=passedFilters(smvcfREF,minCovcutoffREF,maxCovcutoffREF,minMapabilitycutoff,minMQcutoff,minGQcutoff);
	    bool passingFilters=passedFilters(smvcfREF,filtersVCF);
	    //minCovcutoffREF,maxCovcutoffREF,minMapabilitycutoff,minMQcutoff,minGQcutoff);
	    //cout<<rejectFiltersTally()<<endl;
	    //cout<<boolStringify(passingFilters)<<endl;
	    if(!passingFilters){
		continue;
	    }

	    // if(!validOneBP(fieldsEPO[2])){ rejectEPOValidREF++; continue; } //REF in the EPO   
	    // if(!validOneBP(fieldsEPO[4])){ rejectEPOValidANC++; continue; } //chimp/human ancestor in EPO   










		    
	    // if( (smvcfREF.getRef()  != fieldsEPO[2]) ){
	    // 	cerr<<"ERROR: Ref not equal"<<endl;
	    // 	cerr<<smvcfREF<<endl;
	    // 	cerr<<lineFromEPO<<endl;
	    // 	rejectERROR_REF++;
	    // 	continue;
	    // }



		    
	

	    /////////////////////////////////
	    //  PASSED QUALITY FILTERS     //
	    /////////////////////////////////

	    numberOfGoodSites++;

			     
	    if( smvcfREF->isHomozygousREF()  || smvcfREF->isHomozygousALT() ){  countHomo++; }
	    if( smvcfREF->isHeterozygous() ){ countHetero++;	    }

#ifdef DEBUG3
	    cout<<"lineFromREF "<<  *smvcfREF<<endl;
	    // cout<<"lineFromEPO "<<  lineFromEPO<<endl;
#endif


	    //////////////////////////////////////
	    //                                  //
	    //   END  VCF ACTUAL PROCESSING     //
	    //                                  //
	    //////////////////////////////////////





			  
	}else{
#ifdef DEBUG2

	    cout<<"needs line ajusting "<<endl;
#endif


	    //REPOSITING THE FILES AND COUNTER

	    if(coordCurrent < positionREF){ //overshot , repositioning there
		coordCurrent = positionREF;
		continue;
	    }

	    if(coordCurrent == positionREF){ //fine
	    }

	    if(coordCurrent > positionREF){ //running behind
		lineLeftREF=vcfREF.hasData();
		if(lineLeftREF){
		    smvcfREF=vcfREF.getData();
		}
	    }
			    

	    // if(coordCurrent < positionEPO){ //overshot , repositioning there
	    // 	coordCurrent = positionEPO;
	    // 	continue;
	    // }

	    // if(coordCurrent == positionEPO){ //fine
	    // }

	    // if(coordCurrent > positionEPO){ //running behind
	    // 	lineLeftEPO=(rtEPO.readLine( lineFromEPO ));

	    // }
	}
			   
			
		
    }//end of data processing for a given random chromosome (closing while(stayInLoop))

#ifdef DEBUG



    // cout<<"rejectEPOValidREF     "<<rejectEPOValidREF<<endl;
    // cout<<"rejectEPOValidANC     "<<rejectEPOValidANC<<endl;
    cout<<"rejectERROR_REF       "<<rejectERROR_REF<<endl;
    cout<<"rejectREF_genotypeUnknown "<<rejectREF_genotypeUnknown<<endl;

	    
#endif
				
    //we found a satisfactory # of sites
    if(numberOfGoodSites >= minimumNumberOfGoodSites){
	// cout<<chrName<<":"<<startChrCoord<<"-"<<endChrCoord<<"\t"
	//     <<ancestorAllele<<"\t"
	//     <<otherAllele<<"\t"
	//     <<double(ancestorAllele)/double(ancestorAllele+otherAllele)<<endl;
	cout<<chrName<<":"<<startChrCoord<<"-"<<endChrCoord<<"\t"
	    <<countHetero<<"\t"
	    <<countHomo<<"\t"
	    <<double(countHetero)/double(countHetero+countHomo)<<endl;
	return 0;
    }else{
	cerr<<"Failed enough sites good sites for "<<chrName<<":"<<startChrCoord<<"-"<<endChrCoord<<"\t"<<endl;
	cerr<<"Needed "<<minimumNumberOfGoodSites<<" found "<<numberOfGoodSites<<endl;
	return 1;
    }

	    
 }//end of subroutine
