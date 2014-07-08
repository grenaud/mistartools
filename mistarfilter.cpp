/*
 * mistarfilter
 * Date: Mar-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>

#include "utils.h"

#include "MistarParser.h"
#include "GenomicRange.h"
#include "mistarOperations.h"

using namespace std;

int main (int argc, char *argv[]) {
    char mode; //d = noundef, s = segsite, n = no sharing, m=sharing, p = pop subset, b = bedfilter, z=strict no sharing

    string usage=string(""+string(argv[0])+"  [mode]"+
			"\nThis program has the following modes:\n\n"+
			//"\t\t"+"pairdivcompute  To compute divergence pairwise for mistar input\n"+
			"\t\t"+"noundef         No undefined sites for populations\n"+
			"\t\t"+"bedfilter       Filter mistar file using sorted bedfile\n"+
			"\t\t"+"segsite         Just retain segregating sites (or trans./transi)\n"+
			"\t\t"+"popsub          Keep a subset of the populations\n"+
			"\t\t"+"removepop       Remove a subset of the populations\n"+
			"\t\t"+"sharing         Retain sites that share alleles between populations\n"+
			"\t\t"+"nosharing       Retain sites that do not share alleles between populations\n"+
			"\t\t"+"znosharing Retain sites that strickly do not share alleles between populations\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    if( string(argv[1]) == "noundef"){
	mode='d';
	usage=string(string(argv[0])+" "+string(argv[1])+" <options> [mistar file]\n\nDescription:\n\tThis will filter out any site where the allele count is nul 0,0 for both reference and alternative\n"+                                                          
		     "\tOptions:\n"+                                                                                                             
		     "\t\t--allowrootun\tAllow the root to be undefined\n");
    }else{
	if( string(argv[1]) == "segsite"){
	    mode='s';
	    usage=string(string(argv[0])+" "+string(argv[1])+" <options> [mistar file]\n\nDescription:\n\tThis will retain sites where the allele count is greater than 0 for either the reference or alternative for at least one individual\n"+
			 "\tOptions:\n"+
			 "\t\t-ts\tKeep transitions only\n"+
			 "\t\t-tv\tKeep transversions only\n");
	}else{
	    if( string(argv[1]) == "nosharing"){
		mode='n';
		usage=string(string(argv[0])+" "+string(argv[1])+" [mistar file] <comma separated group 1> <comma separated group 2>\n\nDescription:\n\tThis will filter sites where individuals in population group 1 do not share at least one allele with individual in population group 2.\n\tIt requires that the allele count for every individual for both groups be non-zero.\n\tA random allele is picked (biased for allele count) for heterozygous position so do not be surprised if you get different outputs every time.\n\tIn other words, the individuals in the first group have to be all reference and the second all alternative or vice-versa.");
	    }else{
		if( string(argv[1]) == "znosharing"){
		    mode='z';
		    usage=string(string(argv[0])+" "+string(argv[1])+" [mistar file] <comma separated group 1> <comma separated group 2>\n\nDescription:\n\tThis will filter sites where individuals in population group 1 strickly do not share any allele with individual in population group 2.\n\tIt requires that the allele count for every individual for both groups be non-zero.\n\tPlease remember that this will exclude any hetezygous sites\n.");
		    
		}else{
		    if( string(argv[1]) == "sharing"){
			mode='m';
			usage=string(string(argv[0])+" "+string(argv[1])+" [mistar file] <comma separated group 1> <comma separated group 2>\n\nDescription:\n\tThis will only retain sites where every individuals in population group 1 share the same allele(s) as every individual in population group 2.\n\tIt requires that the allele count for every individual for both groups be non-zero.\n\tA random allele is picked (biased for allele count) for heterozygous position so do not be surprised if you get different outputs every time.\n\t");
		    }else{
			
			if( string(argv[1]) == "popsub"){
			    mode='p';
			    usage=string(string(argv[0])+" "+string(argv[1])+" [mistar file] <comma separated group to keep>\n\nDescription:\n\tThis will keep only the population specified in the list.\nPlease note that it will set the alternative allele to 'N' if no population has the alternative allele");
			}else{
			    
			    if( string(argv[1]) == "removepop"){
				mode='r';
				usage=string(string(argv[0])+" "+string(argv[1])+" [mistar file] <comma separated group to remove>\n\nDescription:\n\tThis will remove the population specified in the list.\nPlease note that it will set the alternative allele to 'N' if no population has the alternative allele");
			    }else{
				if( string(argv[1]) == "bedfilter"){
				    mode='b';
				    usage=string(string(argv[0])+" "+string(argv[1])+" [mistar file] <sorted bed file>\n\nDescription:\n\tThis will keep only the positions in the bed file");
				}else{			    
				    cerr << "Invalid mode\nUsage:\n"<<usage<<endl;
				    return 1; 
				}
			    }
			}
		    }
		}
	    }
	}
    }


    if(argc == 2 ||
       (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    unsigned int totalRecords=0;
    unsigned int keptRecords=0;

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    cout<<"#MISTAR"<<endl;    	
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;
    switch(mode){

    case 'b': //bed file
	{

	    string bedFileRegions = string(argv[argc-1]);
	    map< string, vector<GenomicRange> * > * bedRegionsToFilter;
	    bedRegionsToFilter = readBEDSortedfile(bedFileRegions);
	    map< string, unsigned int > * coordinateOfVec=new map< string, unsigned int >();
	    for(map<string,vector<GenomicRange> * >::iterator it = bedRegionsToFilter->begin(); 
		it != bedRegionsToFilter->end(); 
		++it) {
		//cout<<it->first<<endl;
		coordinateOfVec->insert( pair<string ,unsigned int>(it->first,0) );
	    }

	    MistarParser mp   (argv[argc-2]);
	    string chrName="-1";
	    unsigned int previousCoordinate=0;
	    vector<GenomicRange> * currentGr=0;
	    unsigned int currentIndex=0;
	    AlleleRecords * dataRow;
	    bool chrFoundInBed=false;

	    cout<<"#MISTARFILTER:bedfilter "<<bedFileRegions<<endl;
	    cout<<""<<mp.getHeader("#\t")<<"\n";
	    cout<<""<<mp.getDefline()<<"\n";
	    
	    while(mp.hasData()){
		dataRow = mp.getData();
		// cout<<"test\t"<<(*dataRow)<<endl;

		if(dataRow->chr != chrName){
		    if(bedRegionsToFilter->find(dataRow->chr) == bedRegionsToFilter->end() ){
			chrFoundInBed=false;
			currentGr    = 0;		       
		    }else{
			chrFoundInBed=true;
			currentGr    = bedRegionsToFilter->at(dataRow->chr) ;
			currentIndex =   coordinateOfVec->at(dataRow->chr) ;
			if(currentIndex!=0){
			    cerr<<"There seems to be a mix of chromosomes in the mistar file, needs to be sorted chr: "<<chrName<<endl;
			    return 1;
			}
		    }
		    previousCoordinate = dataRow->coordinate;
		    chrName            = dataRow->chr;
		}else{
		    if(previousCoordinate >= dataRow->coordinate){
			cerr<<"There seems to be a unsorted coordinate in the mistar file, needs to be sorted coordinate: "<<previousCoordinate<<endl;
			return 1;
		    }
		}

		if(!chrFoundInBed )
		    goto nextmistarrecord;

	

		while(1){
		    if(currentIndex == currentGr->size())//end of bed file
			goto nextmistarrecord;
		    //ignore
		    //        |---------|
		    // *     
		    if( dataRow->coordinate < currentGr->at(currentIndex).getStartCoord() ){//ignore
			// cout<<"case 1"<<endl;
			goto nextmistarrecord;
		    }

		    //print
		    //        |---------|
		    //            *     
		    if(dataRow->coordinate >= currentGr->at(currentIndex).getStartCoord() &&
		       dataRow->coordinate <= currentGr->at(currentIndex).getEndCoord() ){
			// cout<<"case 2"<<endl;
			cout<<(*dataRow)<<endl;
			keptRecords++;				
			goto nextmistarrecord;
		    }

		    //we are running behind in the bed file
		    //        |---------|
		    //                      *     
		    if(  dataRow->coordinate > currentGr->at(currentIndex).getEndCoord()){
			// cout<<"case 3\t"<<currentIndex<<"\t"<<currentGr->size()<<endl;
			
			if(currentIndex<currentGr->size()){//we move to next iteration			    
			    currentIndex++;
			}else{
			    goto nextmistarrecord; //we have reached the end of the vector, do nothing until next chr
			}
		    }
		}
	    nextmistarrecord:
		// cout<<dataRow->coordinate<<endl;		    

		totalRecords++;
	    }

	    break;	    
	}

    case 'p': //population subset
	{
	    MistarParser mp   (argv[argc-2]);
	    string g1 = string(argv[argc-1]);

	    vector<string> g1v = allTokens(g1,',');
	    // vector<unsigned int>   g1i;
	    vector<bool> flagsPopToAdd = vector<bool>(mp.getPopulationsNames()->size(),false);
	    	    

	    for(unsigned k=0;k<g1v.size();k++){
		bool found=false;
		for(unsigned i=0;i<(mp.getPopulationsNames()->size());i++){
		    if(mp.getPopulationsNames()->at(i) == g1v[k]){
			//g1i.push_back(i);
			found=true;
			//s1i.insert(i);
			flagsPopToAdd[i] = true;
			break;
		    } 
		}
		
		if(!found){
		    cerr<<"Cannot find population "<<g1v[k]<<endl;
		    return 1;
		}
	    }


	    AlleleRecords * dataRow;
	    cout<<"#MISTARFILTER:popsub "<<g1<<endl;
	    cout<<""<<mp.getHeader("#\t")<<"\n";
	    //cout<<""<<mp.getDefline()<<"\n";
	    vector<string> newdefline;
	    //root and anc
	    newdefline.push_back(mp.getPopulationsNames()->at(0));
	    newdefline.push_back(mp.getPopulationsNames()->at(1));

	    for(unsigned i=2;i<(mp.getPopulationsNames()->size());i++){
		// for(unsigned k=0;k<g1i.size();k++){
		//     if(i==g1i[k]){
		if(flagsPopToAdd[i])
		    newdefline.push_back(mp.getPopulationsNames()->at(i));
		//     }
		// }
	    }

	    cout<<"#chr\tcoord\tREF,ALT\t"<<vectorToString(newdefline,"\t")<<endl;

	    while(mp.hasData()){
		dataRow = mp.getData();
		totalRecords++;
		vector<SingleAllele> toprint;
		//root
		toprint.push_back(dataRow->vectorAlleles->at(0));
		//anc
		toprint.push_back(dataRow->vectorAlleles->at(1));
		bool someoneHasAlt=false;//flag to check if someone has the alternative, otherwise, we will set the alt to N

		for(unsigned j=2;j<dataRow->vectorAlleles->size();j++){
		    		    
		    // for(unsigned k=0;k<g1i.size();k++){
		    // 	if(j==g1i[k]){	//found		 
		    // 	    toprint.push_back(dataRow->vectorAlleles->at(j));
		    // 	    someoneHasAlt=someoneHasAlt || (dataRow->vectorAlleles->at(j).getAltCount() != 0);				
		    // 	}
		    // }

		    if(flagsPopToAdd[j]){	       
			toprint.push_back(dataRow->vectorAlleles->at(j));
			someoneHasAlt=someoneHasAlt || (dataRow->vectorAlleles->at(j).getAltCount() != 0);
		    }
		}
		keptRecords++;

		cout<<dataRow->chr<<"\t"<<dataRow->coordinate<<"\t"<<dataRow->ref<<",";
		if(someoneHasAlt)
		    cout<<dataRow->alt;
		else
		    cout<<'N';
		cout<<"\t"<<vectorToString(toprint,"\t")<<endl;
		
	    }

	    break;	    
	}


    case 'r': //remove population subset
	{

	    MistarParser mp   (argv[argc-2]);
	    string g1 = string(argv[argc-1]);

	    vector<string> g1v = allTokens(g1,',');
	    vector<unsigned int>   g1i;

	    	    

	    for(unsigned k=0;k<g1v.size();k++){
		bool found=false;
		for(unsigned i=0;i<(mp.getPopulationsNames()->size());i++){
		    if(mp.getPopulationsNames()->at(i) == g1v[k]){
			g1i.push_back(i);
			//cout<<i<<endl;
			found=true;
			break;
		    } 
		}
		
		if(!found){
		    cerr<<"Cannot find population "<<g1v[k]<<endl;
		    return 1;
		}
	    }

	    AlleleRecords * dataRow;
	    cout<<"#MISTARFILTER:removepop "<<g1<<endl;
	    cout<<""<<mp.getHeader("#\t")<<"\n";

	    vector<string> newdefline;
	    // //root
	    // newdefline.push_back(mp.getPopulationsNames()->at(0));
	    // //anc
	    // newdefline.push_back(mp.getPopulationsNames()->at(1));

	    for(unsigned i=0;i<(mp.getPopulationsNames()->size());i++){
		bool found=false;
		
		for(unsigned k=0;k<g1i.size();k++){
		    if(i==g1i[k]){
			found=true;
		    }
		}
		
		if(!found)
		    newdefline.push_back(mp.getPopulationsNames()->at(i));
	    }
	    
	    cout<<"#chr\tcoord\tREF,ALT\t"<<vectorToString(newdefline,"\t")<<endl;

	    while(mp.hasData()){
		dataRow = mp.getData();
		totalRecords++;
		vector<SingleAllele> toprint;
		// //root
		// toprint.push_back(dataRow->vectorAlleles->at(0));
		// //anc
		// toprint.push_back(dataRow->vectorAlleles->at(1));
		bool someoneHasAlt=false;//flag to check if someone has the alternative, otherwise, we will set the alt to N

		for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){
		    bool found=false;

		    for(unsigned k=0;k<g1i.size();k++){
			if(j==g1i[k]){
			    found=true;
			    break;
			}
		    }
		    
		    if(!found){
			//skip undefined sites 
			toprint.push_back(dataRow->vectorAlleles->at(j));
			someoneHasAlt=someoneHasAlt || (dataRow->vectorAlleles->at(j).getAltCount() != 0);				
		    }
		    
		}
		keptRecords++;

		cout<<dataRow->chr<<"\t"<<dataRow->coordinate<<"\t"<<dataRow->ref<<",";
		if(someoneHasAlt)
		    cout<<dataRow->alt;
		else
		    cout<<'N';
		cout<<"\t"<<vectorToString(toprint,"\t")<<endl;
		
	    }

	    break;	    
	}


    case 'n': //no sharing
	{
	    MistarParser mp   (argv[argc-3]);
	    string g1 = string(argv[argc-2]);
	    string g2 = string(argv[argc-1]);
	    vector<string> g1v = allTokens(g1,',');
	    vector<string> g2v = allTokens(g2,',');
	    vector<unsigned int>   g1i;
	    vector<unsigned int>   g2i;

	    

	    
	    for(unsigned k=0;k<g1v.size();k++){
		for(unsigned j=0;j<g2v.size();j++){
		    if(g1v[k] ==  g2v[j]){
			cerr<<"Cannot specify the same population ("<<g1v[k]<<") in both groups"<<endl;
			return 1;
		    }
		}
	    }

	    for(unsigned k=0;k<g1v.size();k++){
		bool found=false;
		for(unsigned i=0;i<(mp.getPopulationsNames()->size());i++){
		    if(mp.getPopulationsNames()->at(i) == g1v[k]){
			g1i.push_back(i);
			found=true;
			break;
		    } 
		}
		if(!found){
		    cerr<<"Cannot find population "<<g1v[k]<<endl;
		    return 1;
		}
	    }

	    for(unsigned k=0;k<g2v.size();k++){
		bool found=false;
		for(unsigned i=0;i<mp.getPopulationsNames()->size();i++){
		    if(mp.getPopulationsNames()->at(i) == g2v[k]){
			g2i.push_back(i);
			found=true;
			break;
		    } 
		}
		if(!found){
		    cerr<<"Cannot find population "<<g2v[k]<<endl;
		    return 1;
		}
	    }

	    // cout<<vectorToString(g1i)<<endl;
	    // cout<<vectorToString(g2i)<<endl;
	    // exit(1);

	    AlleleRecords * dataRow;
	    cout<<"#MISTARFILTER:nosharing "<<g1<<"-"<<g2<<endl;
	    cout<<""<<mp.getHeader("#\t")<<"\n";
	    cout<<""<<mp.getDefline()<<"\n";

	    while(mp.hasData()){
		int refCountg1=0;
		int refCountg2=0;
		int altCountg1=0;
		int altCountg2=0;

		dataRow = mp.getData();
		totalRecords++;
		for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){
		    
		    
		    for(unsigned k=0;k<g1i.size();k++){
			if(j==g1i[k]){
			    //skip undefined sites 
			    if(dataRow->vectorAlleles->at(j).getRefCount() == 0 && 
			       dataRow->vectorAlleles->at(j).getAltCount() == 0 )
				goto nextiterationnoshare;			       
			    
			    if(dataRow->vectorAlleles->at(j).getRefCount() != 0 && //hetero, increase one at random
			       dataRow->vectorAlleles->at(j).getAltCount() != 0 ){
				//if the random allele picked is the reference, increase the reference counter
				if(dataRow->ref==sampleRandomRefAltAllele(dataRow->ref,
									  dataRow->alt,
									  dataRow->vectorAlleles->at(j).getRefCount(), 
									  dataRow->vectorAlleles->at(j).getAltCount())){
				    refCountg1+=dataRow->vectorAlleles->at(j).getRefCount();
				}else{//otherwise the alt
				    altCountg1+=dataRow->vectorAlleles->at(j).getAltCount();
				}
			    }else{
				refCountg1+=dataRow->vectorAlleles->at(j).getRefCount();
				altCountg1+=dataRow->vectorAlleles->at(j).getAltCount();
			    }

			}
		    }
		    
		    for(unsigned k=0;k<g2i.size();k++){
			if(j==g2i[k]){
			    //skip undefined sites 
			    if(dataRow->vectorAlleles->at(j).getRefCount() == 0 && 
			       dataRow->vectorAlleles->at(j).getAltCount() == 0 )
				goto nextiterationnoshare;			       

			    if(dataRow->vectorAlleles->at(j).getRefCount() != 0 && //hetero, increase one at random
			       dataRow->vectorAlleles->at(j).getAltCount() != 0 ){
				//if the random allele picked is the reference, increase the reference counter
				if(dataRow->ref==sampleRandomRefAltAllele(dataRow->ref,
									  dataRow->alt,
									  dataRow->vectorAlleles->at(j).getRefCount(), 
									  dataRow->vectorAlleles->at(j).getAltCount())){
				    refCountg2+=dataRow->vectorAlleles->at(j).getRefCount();
				}else{//otherwise the alt
				    altCountg2+=dataRow->vectorAlleles->at(j).getAltCount();
				}
			    }else{
				refCountg2+=dataRow->vectorAlleles->at(j).getRefCount();
				altCountg2+=dataRow->vectorAlleles->at(j).getAltCount();
			    }


			}
		    }
		    

		}
		
		//print data row

		if (refCountg1 != 0 && refCountg2 != 0) // if they share reference
		    goto nextiterationnoshare;

		if (altCountg1 != 0 && altCountg2 != 0) // if they share alternative
		    goto nextiterationnoshare;


		cout<<(*dataRow)<<endl;
		keptRecords++;
	    nextiterationnoshare:	    
		continue;
	    }

	    break;	    
	}










    case 'z': //strick no sharing
	{
	    MistarParser mp   (argv[argc-3]);
	    string g1 = string(argv[argc-2]);
	    string g2 = string(argv[argc-1]);
	    vector<string> g1v = allTokens(g1,',');
	    vector<string> g2v = allTokens(g2,',');
	    vector<unsigned int>   g1i;
	    vector<unsigned int>   g2i;

	    

	    
	    for(unsigned k=0;k<g1v.size();k++){
		for(unsigned j=0;j<g2v.size();j++){
		    if(g1v[k] ==  g2v[j]){
			cerr<<"Cannot specify the same population ("<<g1v[k]<<") in both groups"<<endl;
			return 1;
		    }
		}
	    }

	    for(unsigned k=0;k<g1v.size();k++){
		bool found=false;
		for(unsigned i=0;i<(mp.getPopulationsNames()->size());i++){
		    if(mp.getPopulationsNames()->at(i) == g1v[k]){
			g1i.push_back(i);
			found=true;
			break;
		    } 
		}
		if(!found){
		    cerr<<"Cannot find population "<<g1v[k]<<endl;
		    return 1;
		}
	    }

	    for(unsigned k=0;k<g2v.size();k++){
		bool found=false;
		for(unsigned i=0;i<mp.getPopulationsNames()->size();i++){
		    if(mp.getPopulationsNames()->at(i) == g2v[k]){
			g2i.push_back(i);
			found=true;
			break;
		    } 
		}
		if(!found){
		    cerr<<"Cannot find population "<<g2v[k]<<endl;
		    return 1;
		}
	    }

	    // cout<<vectorToString(g1i)<<endl;
	    // cout<<vectorToString(g2i)<<endl;
	    // exit(1);

	    AlleleRecords * dataRow;
	    cout<<"#MISTARFILTER:znosharing "<<g1<<"-"<<g2<<endl;
	    cout<<""<<mp.getHeader("#\t")<<"\n";
	    cout<<""<<mp.getDefline()<<"\n";

	    while(mp.hasData()){
		int refCountg1=0;
		int refCountg2=0;
		int altCountg1=0;
		int altCountg2=0;

		dataRow = mp.getData();
		totalRecords++;
		for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){
		    
		    
		    for(unsigned k=0;k<g1i.size();k++){
			if(j==g1i[k]){
			    //skip undefined sites 
			    if(dataRow->vectorAlleles->at(j).getRefCount() == 0 && 
			       dataRow->vectorAlleles->at(j).getAltCount() == 0 )
				goto nextiterationstrictnoshare;			       
			    
			    if(dataRow->vectorAlleles->at(j).getRefCount() != 0 && //hetero, increase one at random
			       dataRow->vectorAlleles->at(j).getAltCount() != 0 ){
				goto nextiterationstrictnoshare;			       

			    }else{
				refCountg1+=dataRow->vectorAlleles->at(j).getRefCount();
				altCountg1+=dataRow->vectorAlleles->at(j).getAltCount();
			    }

			}
		    }
		    
		    for(unsigned k=0;k<g2i.size();k++){
			if(j==g2i[k]){
			    //skip undefined sites 
			    if(dataRow->vectorAlleles->at(j).getRefCount() == 0 && 
			       dataRow->vectorAlleles->at(j).getAltCount() == 0 )
				goto nextiterationstrictnoshare;			       

			    if(dataRow->vectorAlleles->at(j).getRefCount() != 0 && //hetero, increase one at random
			       dataRow->vectorAlleles->at(j).getAltCount() != 0 ){
				goto nextiterationstrictnoshare;			       
			    }else{
				refCountg2+=dataRow->vectorAlleles->at(j).getRefCount();
				altCountg2+=dataRow->vectorAlleles->at(j).getAltCount();
			    }


			}
		    }
		    

		}
		
		//print data row

		if (refCountg1 != 0 && refCountg2 != 0) // if they share reference
		    goto nextiterationstrictnoshare;

		if (altCountg1 != 0 && altCountg2 != 0) // if they share alternative
		    goto nextiterationstrictnoshare;


		cout<<(*dataRow)<<endl;
		keptRecords++;
	    nextiterationstrictnoshare:	    
		continue;
	    }

	    break;	    
	}



    case 'm': //sharing
	{
	    MistarParser mp   (argv[argc-3]);
	    string g1 = string(argv[argc-2]);
	    string g2 = string(argv[argc-1]);
	    vector<string> g1v = allTokens(g1,',');
	    vector<string> g2v = allTokens(g2,',');
	    vector<unsigned int>   g1i;
	    vector<unsigned int>   g2i;

	    

	    
	    for(unsigned k=0;k<g1v.size();k++){
		for(unsigned j=0;j<g2v.size();j++){
		    if(g1v[k] ==  g2v[j]){
			cerr<<"Cannot specify the same population ("<<g1v[k]<<") in both groups"<<endl;
			return 1;
		    }
		}
	    }	    

	    for(unsigned k=0;k<g1v.size();k++){
		bool found=false;
		for(unsigned i=0;i<(mp.getPopulationsNames()->size());i++){
		    if(mp.getPopulationsNames()->at(i) == g1v[k]){
			g1i.push_back(i);
			found=true;
			break;
		    } 
		}
		if(!found){
		    cerr<<"Cannot find population "<<g1v[k]<<endl;
		    return 1;
		}
	    }

	    for(unsigned k=0;k<g2v.size();k++){
		bool found=false;
		for(unsigned i=0;i<mp.getPopulationsNames()->size();i++){
		    if(mp.getPopulationsNames()->at(i) == g2v[k]){
			g2i.push_back(i);
			found=true;
			break;
		    } 
		}
		if(!found){
		    cerr<<"Cannot find population "<<g2v[k]<<endl;
		    return 1;
		}
	    }


	    AlleleRecords * dataRow;
	    cout<<"#MISTARFILTER:sharing "<<g1<<"-"<<g2<<endl;
	    cout<<""<<mp.getHeader("#\t")<<"\n";
	    cout<<""<<mp.getDefline()<<"\n";

	    while(mp.hasData()){
		int refCountg1=0;
		int refCountg2=0;
		int altCountg1=0;
		int altCountg2=0;

		dataRow = mp.getData();
		totalRecords++;
		for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){
		    
		    
		    for(unsigned k=0;k<g1i.size();k++){
			if(j==g1i[k]){
			    //skip undefined sites 
			    if(dataRow->vectorAlleles->at(j).getRefCount() == 0 && 
			       dataRow->vectorAlleles->at(j).getAltCount() == 0 )
				goto nextiterationshare;			       

			    if(dataRow->vectorAlleles->at(j).getRefCount() != 0 && //hetero, increase one at random
			       dataRow->vectorAlleles->at(j).getAltCount() != 0 ){
				//if the random allele picked is the reference, increase the reference counter
				if(dataRow->ref==sampleRandomRefAltAllele(dataRow->ref,
									  dataRow->alt,
									  dataRow->vectorAlleles->at(j).getRefCount(), 
									  dataRow->vectorAlleles->at(j).getAltCount())){
				    refCountg1+=dataRow->vectorAlleles->at(j).getRefCount();
				}else{//otherwise the alt
				    altCountg1+=dataRow->vectorAlleles->at(j).getAltCount();
				}
			    }else{
				refCountg1+=dataRow->vectorAlleles->at(j).getRefCount();
				altCountg1+=dataRow->vectorAlleles->at(j).getAltCount();
			    }

			}
		    }
		    
		    for(unsigned k=0;k<g2i.size();k++){
			if(j==g2i[k]){
			    //skip undefined sites 
			    if(dataRow->vectorAlleles->at(j).getRefCount() == 0 && 
			       dataRow->vectorAlleles->at(j).getAltCount() == 0 )
				goto nextiterationshare;			       


			    if(dataRow->vectorAlleles->at(j).getRefCount() != 0 && //hetero, increase one at random
			       dataRow->vectorAlleles->at(j).getAltCount() != 0 ){
				//if the random allele picked is the reference, increase the reference counter
				if(dataRow->ref==sampleRandomRefAltAllele(dataRow->ref,
									  dataRow->alt,
									  dataRow->vectorAlleles->at(j).getRefCount(), 
									  dataRow->vectorAlleles->at(j).getAltCount())){
				    refCountg2+=dataRow->vectorAlleles->at(j).getRefCount();
				}else{//otherwise the alt
				    altCountg2+=dataRow->vectorAlleles->at(j).getAltCount();
				}
			    }else{
				refCountg2+=dataRow->vectorAlleles->at(j).getRefCount();
				altCountg2+=dataRow->vectorAlleles->at(j).getAltCount();
			    }
			}
		    }
		    

		}
		
		//print data row

		if (refCountg1 != 0 ){ //if first group has ref
		    if(refCountg2 == 0)  // we must observe ref in the second
			goto nextiterationshare;
		}

		if (refCountg2 != 0 ){ //if second group has ref
		    if(refCountg1 == 0)  // we must observe ref in the first
			goto nextiterationshare;
		}

		if (altCountg1 != 0 ){ //if first group has alt
		    if(altCountg2 == 0) // we must observe alt in the second
			goto nextiterationshare;
		}

		if (altCountg2 != 0 ){ //if second group has alt
		    if(altCountg1 == 0) // we must observe alt in the first
			goto nextiterationshare;
		}


		cout<<(*dataRow)<<endl;
		keptRecords++;
	    nextiterationshare:	    
		continue;
	    }

	    break;	    
	}






    case 's': //just segregating sites
	{

	    bool onlyTS=false;
	    bool onlyTV=false;


	    for(int i=0;i<(argc-1);i++){ 
		if(string(argv[i]) == "-ts"){
		    onlyTS=true;
		}

		if(string(argv[i]) == "-tv"){
		    onlyTV=true;
		}

	    }
	    
	    if(onlyTS == true &&
	       onlyTV == true ){
		cerr<<"Cannot select on transitions and transversions at the same time"<<endl;
		return 1;
	    }
	       


	    MistarParser mp (argv[argc-1]);
	    AlleleRecords * dataRow;


	    cout<<"#MISTARFILTER:segsite"<<endl;
	    cout<<""<<mp.getHeader("#\t")<<"\n";
	    cout<<""<<mp.getDefline()<<"\n";

	    while(mp.hasData()){
		dataRow = mp.getData();
		totalRecords++;
		bool haveSeenRef=false;
		bool haveSeenAlt=false;

		for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){
		    if( (dataRow->vectorAlleles->at(j).getRefCount() != 0) ){
			haveSeenRef=true;
		    }
		    if( (dataRow->vectorAlleles->at(j).getAltCount() != 0) ){
			haveSeenAlt=true;
		    }
		}
		//print data row if have seen at least once both
		if(haveSeenRef && haveSeenAlt){
		    bool isTrans = isPotentialTransition(dataRow->ref,dataRow->alt);
		    if(onlyTV || onlyTS){

			if(onlyTV && !isTrans){
			    cout<<(*dataRow)<<endl;
			}
			
			if(onlyTS && isTrans){
			    cout<<(*dataRow)<<endl;
			}


		    }else{		    
			cout<<(*dataRow)<<endl;
		    }
		    keptRecords++;
		}
	    }

	    break;	    
	}


    case 'd': //no undefined (just zeros) site 
	{
	    unsigned int firstPopInd=0;
	    string mistarFileToUn;
	    if(string(argv[argc-2]) == "--allowrootun"){
		firstPopInd=2;		
	    }

	    mistarFileToUn = string(argv[argc-1]);
	    MistarParser mp (mistarFileToUn);
	    AlleleRecords * dataRow;


	    cout<<"#MISTARFILTER:noundef"<<endl;
	    cout<<""<<mp.getHeader("#\t")<<"\n";
	    cout<<""<<mp.getDefline()<<"\n";

	    while(mp.hasData()){
		dataRow = mp.getData();
		totalRecords++;

		for(unsigned j=firstPopInd;j<dataRow->vectorAlleles->size();j++){
		    //undefined site
		    if( (dataRow->vectorAlleles->at(j).getRefCount() == 0) &&
			(dataRow->vectorAlleles->at(j).getAltCount() == 0) ){
			goto nextiterationnoundef;
		    }
		}

		cout<<(*dataRow)<<endl;
		keptRecords++;

	    nextiterationnoundef:	    
		continue;
	    }

	    break;
	}



    default:
	cerr << "Unknown mode "<<mode<<endl;
	return 1;       
	break;	    

    }
    cerr<<"Program "<<argv[0]<<" wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully";

    return 0;
}

