/*
 * ms2mistar
 * Date: May-14-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <inttypes.h>
#include <sys/time.h> //for srand

#include "MSParser.h"
#include "MSobject.h"
#include "SingleAllele.h"

#include "utils.h"

using namespace std;


//TODO 
// realistic base pair #
// ability to pair populations


char randomBase(){
    return "ACGT"[ rand()%4 ];
}

char randomBaseNotThat(char c){
    while(1){
	char toreturn="ACGT"[ rand()%4 ];
	if(toreturn!=c)
	    return toreturn;
    }
}



int main (int argc, char *argv[]) {

    timeval time;
    gettimeofday(&time, NULL);
    srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );
    srandCalled=true;

    int numberOfIndividuals=-1;
    //define options : ms file, which individuals are "root"
    //string usage=string("\t"+string(argv[0])+"  <mistar file> <correspondance individuals to pop> <size of chromosome>\n"+
    string usage=string("\t"+string(argv[0])+"  [mistar file] [correspondance individuals to pop] [size of chromosome]\n"+
			"This program converts ms output into a mistar matrix\nThe correspondence has to have the following format\npop1:individual1,individual2-npop2:individual3,individual4\n"+
			"The size of the chromsome is the parameter used as -r\n"+
			"\nExample:\n\t"+
			//string(argv[0])+" input.msout root:1,2-pop2:3-4 10000\n"+
			string(argv[0])+" input.msout root:1,2-pop2:3-4 10000\n"+
			//"\nThe genome length is the value you passed to -r"+
			"\nThe first individual is individual #1 (not zero)");

     
    if(argc != 4 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    cout<<"#MISTAR"<<endl;
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;

    cout<<"#chr\tcoord\tREF,ALT\troot\tanc\t";
    

    MSParser msp                                 ( string(argv[argc-3]) );
    string individuals                           = string(argv[argc-2]);
    unsigned int sizeChr=destringify<unsigned int>(string(argv[argc-1]));
    cerr<<sizeChr<<endl;

    //Parsing the populations argument
    vector<string> popNames;
    bool rootWasSpecified=false;
    vector<string> tempind=allTokens(individuals,'-');
    int maxInd=0;
    int popFound=0;
    vector<string> indString;

    //counting # of populations 
    for(unsigned int i=0;i<tempind.size();i++){
	popFound++;

	//cerr<<tempind[i]<<endl;
	vector<string> tempind2   = allTokens( tempind[i],':');
	indString.push_back(tempind2[0]);

	if(tempind2[0]=="root" ){
	    rootWasSpecified=true;
	    if(i != 0){
		cerr<<"If you specify the root, it must be the first pop"<<endl;
		return 1;
	    }	
	}
	popNames.push_back(tempind2[0]);
	cerr<<"pop: "<<tempind2[0]<<"\t";
	//which individual
	vector<string> popIndices = allTokens(tempind2[1],',');
	for(unsigned int j=0;j<popIndices.size();j++){
	    int x=destringify<int>(popIndices[j]);
	    cerr<<" "<<x<<" ";
	    if( maxInd < x  ){
		maxInd=x;
	    }
	}
	cerr<<endl;
    }

    vector<unsigned int> ind2pop (maxInd,0);
    unsigned int currentPopIndex=0;

    for(unsigned int i=0;i<tempind.size();i++){
	vector<string> tempind2   = allTokens( tempind[i],':');
	vector<string> popIndices = allTokens(tempind2[1],',');
	for(unsigned int j=0;j<popIndices.size();j++){
	    int x=destringify<int>(popIndices[j])-1;
	    ind2pop[ x ]=currentPopIndex;
	}
	currentPopIndex++;
    }

    // cout<<vectorToString(ind2pop)<<endl;
    // return 1;
    //cout<<string(argv[argc-1])<<endl;

    

    vector<unsigned int> rootInd;
   

    unsigned int chromosome=1;
    // uint64_t  coordinate;
    bool firstRecord=true;
    uint64_t records=0;
    // vector<bool> isRoot;
    // cout<<"ok"<<endl;

    for(int numRec=0;numRec<msp.numberOfRecords();numRec++){
	const MSobject * mso = msp.getMSObj(numRec);

	if(firstRecord){
	    numberOfIndividuals=mso->getNumberIndividuals();	    

	    if(numberOfIndividuals != maxInd){
		cerr<<"The number of individuals you entered "<<maxInd<<" is different from the one we found "<<numberOfIndividuals<<endl;
		return 1;
	    }
	    vector<string> toprintPop=indString;
	    if(rootWasSpecified){
		toprintPop.erase(toprintPop.begin());//the root should be the first element
	    }
	    cout<<vectorToString(toprintPop,"\t")<<endl;
	    // return 1;
	    firstRecord=false;
	}

	if(numberOfIndividuals != mso->getNumberIndividuals()){
	    cerr<<"Uneven number of individuals"<<endl;
	    return 1;
	}
	
	vector<const bool *> indAllele;
	for(int i=1;i<=mso->getNumberIndividuals();i++){
	    indAllele.push_back(mso->getBoolArrayIndividual(i));
	}

	unsigned int lastCoord=0;

	for(unsigned int k=0;k<mso->getNumberSegSites();k++){

	    //printing non-seg sites as homo ref sites
	    for(unsigned int tempCoord=lastCoord;
		tempCoord<( (mso->getPositions()->at(k)));
		tempCoord++){
		// cout<<tempCoord<<endl;
		cout<<chromosome<<"\t"<< (tempCoord)<<"\t";
		 char randTempBase=randomBase();
		 cout<<randTempBase<<",N"<<"\t";
		 vector<SingleAllele> saToprint (popFound+1,SingleAllele());

		 for(int j=0;j<(popFound+1);j++){
		    saToprint[j].addRefCount(1);
		 }
		 cout<<vectorToString(saToprint,"\t")<<endl;
	    }
		
	    
	    cout<<chromosome<<"\t"<< (mso->getPositions()->at(k))<<"\t";
	    //cerr<<(mso->getPositions()->at(k))<<endl;
	    lastCoord=(mso->getPositions()->at(k));
	    char ref=randomBase();
	    char alt=randomBaseNotThat(ref);
	    cout<<ref<<","<<alt<<"\t";
	    //cout<<endl;
	    
	    SingleAllele rootAl;
	    SingleAllele ancAl;
	    ancAl.addRefCount(1);
			
	    vector<SingleAllele> saToprint (popFound,SingleAllele());

	    for(int j=0;j<numberOfIndividuals;j++){
		unsigned int indexOfPop=ind2pop[j];
		// //if(isRoot[j]){
		if(indAllele[j][k] == 0)//anc
		    saToprint[indexOfPop].addRefCount(1);
		else
		    saToprint[indexOfPop].addAltCount(1);
	    }
	    
	    if(rootWasSpecified){
		cout<<saToprint[0]<<"\t"<<ancAl<<"\t";
		saToprint.erase(saToprint.begin());
		cout<<vectorToString(saToprint,"\t")<<endl;
	    }else{
		cout<<rootAl<<"\t"<<ancAl<<"\t"<<vectorToString(saToprint,"\t")<<endl;
	    }
	    records++;
	}//for each seg site

	// for(int i=0;i<numberOfIndividuals;i++){
	//     namesToUse.push_back("i"+stringify(i));
	// }

	//print empty records till end
	for(unsigned int tempCoord=lastCoord;
	    tempCoord<sizeChr;
	    tempCoord++){
	    //cerr<<tempCoord<<endl;
	    cout<<chromosome<<"\t"<< (tempCoord)<<"\t";
	    char randTempBase=randomBase();
	    cout<<randTempBase<<",N"<<"\t";
	    vector<SingleAllele> saToprint (popFound+1,SingleAllele());
	    
	    for(int j=0;j<(popFound+1);j++){
		saToprint[j].addRefCount(1);
	    }
	    cout<<vectorToString(saToprint,"\t")<<endl;
	}


	chromosome++;
    }

    cerr<<"Program terminated successfully, wrote "<<records<<" records"<<endl;
    
    return 0;
}

