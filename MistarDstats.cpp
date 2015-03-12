/*
 * 
 * MistarDstats
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "MistarDstats.h"
#include "MistarParser.h"



void dstatsMistar(string filename){
    
   MistarParser mp (filename);

   unsigned int numberOfPopulations=mp.getPopulationsNames()->size()+1;//+1 for the human reference

   //DivergenceResult  divergenceResults[numberOfPopulations][numberOfPopulations];
   //DstatResult dstatResults[numberOfPopulations][numberOfPopulations][numberOfPopulations];
   vector< vector< vector<DstatResult> >  > dstatResults;
   
   dstatResults.resize(numberOfPopulations);
   
   for(unsigned int i=0;i<numberOfPopulations;i++)
       dstatResults[i].resize(numberOfPopulations);
   
   for(unsigned int i=0;i<numberOfPopulations;i++)
       for(unsigned int j=0;j<numberOfPopulations;j++)
	   dstatResults[i][j].resize(numberOfPopulations);

   // struct DstatResultArray{
   //     DstatResult dstatResults[numberOfPopulations][numberOfPopulations][numberOfPopulations];
   // };

   AlleleRecords * currentRow;   
   vector<AlleleRecords> segregatingSites;

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
	   cerr<<"MistarDstats.cpp  dstats() Problem for line "<<currentRow->chr<<" "<<currentRow->coordinate<<" reference = "<<currentRow->ref<<" is not resolved"<<endl;
	   exit(1);
       }
       
       
       //first one is the ancestral
       //last one is the human reference
       char sampledAllele [ numberOfPopulations]; //array of sampled alleles
       bool cpgForPop     [ numberOfPopulations]; //array of flags to say if the current record is cpg

       //initialize the sampledAllele and cpgForPop
       if(!isResolvedDNA(currentRow->alt)){ // if one of A,C,G,T
	   continue; //next iteration, we need a valid alternative allele
       }

       // cout<<"1 "<<currentRow->vectorAlleles->size()<<endl;   
       // cout<<segregatingSites[segregatingSites.size()-1].vectorAlleles->size()<<endl;
       // for(unsigned int m=0;m<segregatingSites.size();m++){
       // 	   cout<<segregatingSites[m]<<endl;
       // }
       // cout<<"2 "<<currentRow->vectorAlleles->size()<<endl;
       //double check, already checked in MistarParser 
       if(currentRow->vectorAlleles->size() != (numberOfPopulations-1)){
	   cerr<<"MistarDstats.cpp  pairwiseDivergence() Problem for line "<<currentRow->chr<<" "<<currentRow->coordinate<<" wrong number of columns"<<endl;
	   exit(1);
       }

       for(unsigned int i=0;i<currentRow->vectorAlleles->size();i++){
	   if(i == 0 ){ //the root can be absent, need to check	      
	       //if the allele count is unknown for one, skip
	       if(currentRow->vectorAlleles->at(i).getRefCount() ==  0 &&
		  currentRow->vectorAlleles->at(i).getAltCount() ==  0){
		   goto SKIPTONEXTITERATION;
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


       //the first is the ancestral
       if(sampledAllele[0] == 'N'){//the root has an unresolved allele, skip
	   goto SKIPTONEXTITERATION;       	   
	   continue;
       }
       // cout<<"record "<<alleleRecordsAsString(*currentRow)<<endl;
       // cout<<"record "<<*currentRow<<endl;


       segregatingSites.push_back(*currentRow);//copy constructor
       // cout<<segregatingSites.size()<<endl;

       
       //for each population, except the ancestral at index 0
       for(unsigned i=1;i<numberOfPopulations;i++){

	   //for each population, except the ancestral at index 0
       	   for(unsigned j=1;j<numberOfPopulations;j++){	       
	       //for each population, except the ancestral at index 0
	       for(unsigned k=1;k<numberOfPopulations;k++){	       
		   

		   //skip when the populations are identical
		   if(i==j || i==k || j==k)
		       continue;
		   // cout<<populationNames->at(j)<<"-"<<populationNames->at(k)<<"@"<< populationNames->at(i)<<endl;
		   
		   bool dstval = computeDstat(sampledAllele[0], //root
					      sampledAllele[i], //derived
					      sampledAllele[j], //ind 1
					      sampledAllele[k], //ind 2
					      (cpgForPop[j] || cpgForPop[k]), //only look at j and k for CpG
					      &(dstatResults[i][j][k]) );
		   
		   // if(dstval && j==2 && k==3 && i==6){
		       
		   //     cout<<currentRow->chr<<"\t"<<currentRow->coordinate<<"\t";
		   //     for(unsigned int q=0;q<numberOfPopulations;q++){
		   // 	   cout<<sampledAllele[q];
		   //     }
		   //     cout<<endl;
		   // }
	       
	       }
	   }
       }

   SKIPTONEXTITERATION:
       continue;
       
   }//while the parser has data

   cerr<<"done with main calculations, doing bootstraps "<<endl;
   //doing bootstraping

   unsigned int numberOfBootstraps=1000;

   // cout<<"done "<<endl;
   timeval time;
   gettimeofday(&time, NULL);
   srand48(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );
   // DstatResult dstatResultsBoot[numberOfBootstraps][numberOfPopulations][numberOfPopulations][numberOfPopulations];
   list<vector< vector< vector<DstatResult> >  > >   dstatResultsList;

   for(unsigned int m=0;m<numberOfBootstraps;m++){
       // cout<<"boot "<<m<<endl;
       //vector<AlleleRecords> currentSet;       
       
       vector< vector< vector<DstatResult> >  > dstatResultsBoot;       
       dstatResultsBoot.resize(numberOfPopulations);  
       for(unsigned int i=0;i<numberOfPopulations;i++)
	   dstatResultsBoot[i].resize(numberOfPopulations);       
       for(unsigned int i=0;i<numberOfPopulations;i++)
	   for(unsigned int j=0;j<numberOfPopulations;j++)
	       dstatResultsBoot[i][j].resize(numberOfPopulations);

       for(unsigned int n=0;n<segregatingSites.size();n++){
	   unsigned int toReturn = (unsigned int)(mrand48()); //should be large enough
	   unsigned int indexrand = toReturn%segregatingSites.size();
	   const AlleleRecords * touse = &segregatingSites[ indexrand ]; //pick random seg site
	   // cout<<"touse "<<indexrand<<"\t"<<*touse<<endl;
	   //no need to do the additional checks on the root and so forth since the push_back occurs 
	   //after the checks in the original loop
	   
	   char sampledAllele [ numberOfPopulations]; //array of sampled alleles
	   bool cpgForPop     [ numberOfPopulations]; //array of flags to say if the current record is cpg

	   for(unsigned int i=0;i<touse->vectorAlleles->size();i++){   
	       //plus one for the human allele in sampledAllele and cpgForPop
	       sampledAllele[i] =  sampleRandomRefAltAllele(touse->ref,touse->alt,
							    touse->vectorAlleles->at(i).getRefCount(),
							    touse->vectorAlleles->at(i).getAltCount());
	       cpgForPop[i] = touse->vectorAlleles->at(i).getIsCpg();
	       //cout<<"sampledAllele["<<i<<"] "<<sampledAllele[i]<<endl;
	   }
	   sampledAllele[numberOfPopulations-1] = touse->ref;
	   cpgForPop[ numberOfPopulations-1]    = touse->vectorAlleles->at(0).getIsCpg();//set the cpg to the ancestral CpG flag


	   // cout<<"---------"<<numberOfPopulations<<endl;


	   //for each population, except the ancestral at index 0
	   for(unsigned int i=1;i<numberOfPopulations;i++){

	       //for each population, except the ancestral at index 0
	       for(unsigned int j=1;j<numberOfPopulations;j++){	       
		   //for each population, except the ancestral at index 0
		   for(unsigned int k=1;k<numberOfPopulations;k++){	       


		       //skip when the populations are identical
		       if(i==j || i==k || j==k)
			   continue;

		       bool dstval = computeDstat(sampledAllele[0], //root
						  sampledAllele[i], //derived
						  sampledAllele[j], //ind 1
						  sampledAllele[k], //ind 2
						  (cpgForPop[j] || cpgForPop[k]), //only look at j and k for CpG
						  &(dstatResultsBoot[i][j][k]) );	  
		       
		      
		   }
	       }
	   }//done for each triplet

       }//done for each seg site
       //store dstatResultsBoot
       
       // for(unsigned i=1;i<numberOfPopulations;i++){
       //     for(unsigned j=1;j<numberOfPopulations;j++){	    
       // 	   for(unsigned k=1;k<numberOfPopulations;k++){	       

       // 	       //for each population, except the ancestral at index 0

       // 	   //skip when the allele is identical
       // 	   if(i==j || i==k || j==k)
       // 	       continue;	   
       // 	   cout<<populationNames->at(j)<<"-"<<populationNames->at(k)<<"@"<< populationNames->at(i)  <<"\t"<<dstatResultsBoot[i][j][k]<<endl;
       // 	   }
       //     }
       // }
       // break;
       dstatResultsList.push_back(dstatResultsBoot);
   }//done for each bootstrap





   //printing      
   for(unsigned int  i=1;i<numberOfPopulations;i++){
       for(unsigned int j=1;j<numberOfPopulations;j++){	    
	   for(unsigned int k=1;k<numberOfPopulations;k++){	       	       
   	       //for each population, except the ancestral at index 0
	       //skip when the allele is identical
	       if(i==j || i==k || j==k)
		   continue;	   

	       cout<<populationNames->at(j)<<"-"<<populationNames->at(k)<<"@"<< populationNames->at(i)  <<"\t"<<dstatResults[i][j][k].printWithBootstrap(dstatResultsList,i,j,k,numberOfBootstraps)<<endl;
	   }
       }
   }

   delete(populationNames);
}

