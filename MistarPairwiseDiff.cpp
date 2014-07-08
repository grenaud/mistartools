/*
 * 
MistarPairwiseDiff
 * Date: Jan-26-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "MistarPairwiseDiff.h"



AllPairDistanceResult *  pairwiseDifferences_calc( MistarParser * mp, string chrNameRepo,unsigned int startCoordRepo,unsigned int endCoordRepo  ){
    unsigned int numberOfPopulations=mp->getPopulationsNames()->size();
    
    AllPairDistanceResult * adr =new AllPairDistanceResult(numberOfPopulations,*mp->getPopulationsNames());
    
    AlleleRecords * currentRow;
    //const vector<string> * populationNames = mp.getPopulationsNames();

    cerr<<"NOTE: Please rememeber that the distance to the root is independent of the distance to the remaining populations, in other words, if the root information is not available, the distance between the remaining individuals will still get printed"<<endl;

   //reposition tabix iterator if the calling function specified regions
   if(chrNameRepo       !=  "NONE" && 
      startCoordRepo    !=  0      && 
      endCoordRepo      !=  0      ){
       mp->repositionIterator(chrNameRepo,startCoordRepo,endCoordRepo);
   }

      


   while(mp->hasData()){
       currentRow = mp->getData();
       // cout<<*currentRow<<endl;
       if(!isResolvedDNA(currentRow->ref) ){
	   cerr<<"MistarPairwiseDiff.cpp  pairwiseDifferences() Problem for line "<<currentRow->chr<<" "<<currentRow->coordinate<<" reference = "<<currentRow->ref<<" is not resolved"<<endl;
	   exit(1);
       }


       char sampledAllele[ numberOfPopulations]; //array of sampled alleles
       bool cpgForPop    [ numberOfPopulations]; //array of flags to say if the current record is cpg


       //initialize the sampledAllele and cpgForPop
       if(isResolvedDNA(currentRow->alt)){ // if one of A,C,G,T
	   for(unsigned i=0;i<currentRow->vectorAlleles->size();i++){
	       if(i == 0 || i == 1 ){ //the root can be absent
		   if(currentRow->vectorAlleles->at(i).getRefCount() ==  0 &&
		      currentRow->vectorAlleles->at(i).getAltCount() ==  0){
		       sampledAllele[i] = 'N'; //dummy value
		       cpgForPop[i]     = false;
		       continue;
		   }	       
	       }
	       sampledAllele[i] =  sampleRandomRefAltAllele(currentRow->ref,
							    currentRow->alt,
							    currentRow->vectorAlleles->at(i).getRefCount(),
							    currentRow->vectorAlleles->at(i).getAltCount());
	       cpgForPop[i] = currentRow->vectorAlleles->at(i).getIsCpg();
	   }
	   
       }else{ //if the alt is 'N'


	   if( currentRow->alt != 'N' ){ 
	       cerr<<"MistarPairwiseDiff.cpp  pairwiseDifferences() Problem for line "<<currentRow->chr<<" "<<currentRow->coordinate<<" alternative = "<<currentRow->alt<<" is not resolved"<<endl;
	       exit(1);
	   }

	    for(unsigned i=0;i<currentRow->vectorAlleles->size();i++){
	       if(i == 0 || i == 1){ //the root or ancestor can be absent
		   if(currentRow->vectorAlleles->at(i).getRefCount() ==  0 &&
		      currentRow->vectorAlleles->at(i).getAltCount() ==  0){
		       sampledAllele[i] = 'N'; //dummy value
		       cpgForPop[i]     = false;
		       continue;
		   }	       
	       }

	       //they are all identical
	       sampledAllele[i] =  currentRow->ref;
	       cpgForPop[i] = currentRow->vectorAlleles->at(i).getIsCpg();
	   }
       }



       //for each population
       for(unsigned i=0;i<currentRow->vectorAlleles->size();i++){

	   //the first is the chimp, second is the ancestral
	   if(i == 0 || i==1 ){
	       if(sampledAllele[i] == 'N'){//the root has an unresolved allele, skip
		   continue;
	       }
	   }

       	   for(unsigned j=(i+1);j<currentRow->vectorAlleles->size();j++){	       

	       int allePairIndex=allelePair2Int(sampledAllele[i],sampledAllele[j]);

	       adr->distResults->at(i)->at(j).all.addAllelePair(allePairIndex);

	       //is cpg ?
	       if(cpgForPop[i] || cpgForPop[j])
		   adr->distResults->at(i)->at(j).onlyCpg.addAllelePair(allePairIndex);
	       else
		   adr->distResults->at(i)->at(j).noCpg.addAllelePair(allePairIndex);

	       //it's a mutation, transversion or transition?
	       if(sampledAllele[i] != sampledAllele[j]){
		   if(isPotentialTransition(sampledAllele[i],sampledAllele[j]))
		       adr->distResults->at(i)->at(j).transitions.addAllelePair(allePairIndex);
		   else
		       adr->distResults->at(i)->at(j).transversions.addAllelePair(allePairIndex);
	       }
	       
	   }//each pop 2
       }//end each pop 1

	       
   }

   // vector< vector<DistanceResult> > vecDist;
   // vector<string> namesToUse;
   // for(unsigned i=0;i<numberOfPopulations;i++){
   //     if(i==1) //no anc
   // 	   continue;
   //     namesToUse.push_back(mp.getPopulationsNames()->at(i));

   //     vector<DistanceResult> temp;
   //     for(unsigned j=(0);j<numberOfPopulations;j++){	       
   // 	   if(j==1) //no anc
   // 	       continue;

   // 	   if(j>=(i+1)){
   // 	       if(!computeTree){
   // 		   cout<<populationNames->at(i)<<"-"<<populationNames->at(j)<<endl;
   // 		   cout<<adr->distResults->at(i)->at(j)<<endl;
   // 	       }
   // 	       temp.push_back(adr->distResults->at(i)->at(j));
   // 	   }else{
   // 	       DistanceResult empty;
   // 	       temp.push_back(empty);
   // 	   }
   //     }
   //     vecDist.push_back(temp);
   // }

   // cerr<<"calling"<<endl;

   return adr;
}


AllPairDistanceResult * pairwiseDifferences(string filename){
    
   MistarParser mp (filename);

   return pairwiseDifferences_calc(&mp);
}


