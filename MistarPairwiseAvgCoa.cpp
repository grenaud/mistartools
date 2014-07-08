/*
 * 
 * MistarPairwiseAvgCoa
 * Date: Feb-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "MistarPairwiseAvgCoa.h"
#include "MistarParser.h"



void pairwiseAvgCoa(string filename,bool allowUndefined){
    
   MistarParser mp (filename);

   unsigned int numberOfPopulations=mp.getPopulationsNames()->size()+1;//+1 for the human reference

   AvgCoaResult  divergenceResults[numberOfPopulations][numberOfPopulations];
   AlleleRecords * currentRow;
   vector<string> * populationNames = new vector<string>();
   for(unsigned i=0;i<mp.getPopulationsNames()->size();i++){
       populationNames->push_back(  mp.getPopulationsNames()->at(i) );
   }
   populationNames->push_back("href");
   // for(unsigned i=0;i<populationNames->size();i++){
   //     cout<<i<<"\t"<<populationNames->at(i)<<endl;
   // }
   // exit(1);

   while(mp.hasData()){
       currentRow = mp.getData();
       if(!isResolvedDNA(currentRow->ref) ){
	   cerr<<"MistarPairwiseAvgCoa.cpp  pairwiseAvgCoa() Problem for line "<<currentRow->chr<<" "<<currentRow->coordinate<<" reference = "<<currentRow->ref<<" is not resolved"<<endl;
	   exit(1);
       }

       //first one is the ancestral
       //last one is the human reference
       char sampledAllele[ numberOfPopulations]; //array of sampled alleles
       bool cpgForPop    [ numberOfPopulations]; //array of flags to say if the current record is cpg

       //initialize the sampledAllele and cpgForPop
       if(!isResolvedDNA(currentRow->alt)){ // if one of A,C,G,T
       	   continue; //next iteration, we need a valid alternative allele otherwise there is no impact on divergence calculation, this speeds it up
       }

       //double check, already checked in MistarParser 
       if(currentRow->vectorAlleles->size() != (numberOfPopulations-1)){
	   cerr<<"MistarPairwiseAvgCoa.cpp  pairwiseAvgCoa() Problem for line "<<currentRow->chr<<" "<<currentRow->coordinate<<" wrong number of columns"<<endl;
	   exit(1);
       }

       // cout<<"state 2"<<endl;
       // cout<<currentRow->coordinate<<endl;

       //start at 1 for ancestral
       for(unsigned i=1;i<currentRow->vectorAlleles->size();i++){
	   if(i == 1 ){ //the root can be absent, need to check	      
	       //if the allele count is unknown for both, skip
	       if(currentRow->vectorAlleles->at(i).getRefCount() ==  0 &&
		  currentRow->vectorAlleles->at(i).getAltCount() ==  0){
		   goto SKIPTONEXTITERATION;
	       }
	   } 

	   if(allowUndefined){
	       if(currentRow->vectorAlleles->at(i).getRefCount() ==  0 &&
		  currentRow->vectorAlleles->at(i).getAltCount() ==  0){
		   sampledAllele[i] = 'N';
		   cpgForPop[i]     = false;
		   continue;
	       }
	   }

	   //plus one for the human allele in sampledAllele and cpgForPop
	   sampledAllele[i] =  sampleRandomRefAltAllele(currentRow->ref,currentRow->alt,
							currentRow->vectorAlleles->at(i).getRefCount(),
							currentRow->vectorAlleles->at(i).getAltCount());
	   cpgForPop[i] = currentRow->vectorAlleles->at(i).getIsCpg();
       }

       //storing the human refernce
       sampledAllele[numberOfPopulations-1] = currentRow->ref;
       cpgForPop[ numberOfPopulations-1]    = currentRow->vectorAlleles->at(0).getIsCpg();//set the cpg to the ancestral CpG flag

       //debug
       // cout<<"state 2"<<endl;
       // cout<<currentRow->coordinate<<endl;
       // for(unsigned int ind=0;ind<numberOfPopulations;ind++)
       // 	   cout<<"sampledAllele["<<ind<<"]\t"<<sampledAllele[ind]<<endl;
       // //exit(1);


       //the first is the ancestral, we should never reach that state given the check above
       if(sampledAllele[1] == 'N'){//the ancestral has an unresolved allele, skip
	   goto SKIPTONEXTITERATION;
	   continue;
       }
       // cout<<"record "<<alleleRecordsAsString(*currentRow)<<endl;


       //for each population, except the root/ancestral at index 0,1
       for(unsigned i=2;i<numberOfPopulations;i++){

	   //for each population, except the root/ancestral at index 0,1
       	   for(unsigned j=2;j<numberOfPopulations;j++){	       
	       //skip when the allele is identical
	       if(i==j)
		   continue;
	       	      

	       if(allowUndefined){//if one has undefined allele
		   if(sampledAllele[i] == 'N')
		       continue;
		   if(sampledAllele[j] == 'N')
		       continue;
	       }

	       computeDiv(sampledAllele[1], //0 is root, 1 is ancestral
			  sampledAllele[i],
			  sampledAllele[j],
			  (cpgForPop[i] || cpgForPop[j]),
			  &(divergenceResults[i][j])  );
	       
       	   }
       }

   SKIPTONEXTITERATION:
       continue;
       
   }//while the parser has data




   //printing      
   for(unsigned i=2;i<numberOfPopulations;i++){
       for(unsigned j=2;j<numberOfPopulations;j++){	       
	   //skip when the allele is identical
	   if(i==j)
	       continue;	   
	   cout<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"\t"<<divergenceResults[i][j]<<endl;
       }
   }

   delete(populationNames);
}

