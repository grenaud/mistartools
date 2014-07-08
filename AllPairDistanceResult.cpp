/*
 * AllPairDistanceResult
 * Date: May-01-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "AllPairDistanceResult.h"

AllPairDistanceResult::AllPairDistanceResult(int numberTaxon,vector<string> _populationNames){
    //distResults = new   vector< vector<DistanceResult> * > (numberTaxon, new vector<DistanceResult> (numberTaxon));
    distResults = new   vector< vector<DistanceResult> * > ();

    for(unsigned int i=0;i<_populationNames.size();i++){
	if(_populationNames[i].find("-") != string::npos ){
	    cerr<<"ERROR: AllPairDistanceResult: Population name cannot contain a dash "<<_populationNames[i]<<endl;
	    exit(1);
	}
    }

    for(int i=0;i<numberTaxon;i++){
	vector<DistanceResult> * tempdist = new  vector<DistanceResult>  ();
	for(int j=0;j<numberTaxon;j++){
	    DistanceResult empty;
	    tempdist->push_back(empty);
	}
	distResults->push_back(tempdist);
    }

    numberOfPopulations=numberTaxon;
    populationNames=_populationNames;//copy
}


AllPairDistanceResult::~AllPairDistanceResult(){
    delete distResults;
}

AllPairDistanceResult::AllPairDistanceResult(const string stringForm){
    // ifstream  finSt(stringForm.c_str());
    // string    line;
    
    // while(getline(finSt, line)){
    //   //current line of text is in file_line, not including the \n 
    // 	cout<<"line "<<line<<endl;
    // }
    vector<string> allLines=allTokens(stringForm,'\n');
    distResults = new   vector< vector<DistanceResult> * > ();

    numberOfPopulations=0;
    bool previousEmpty=true;
    populationNames.push_back("root");
    populationNames.push_back("anc");

    for(unsigned int i=0;i<allLines.size();i++){
	// cout<<allLines[i].empty()<<"\t"<<allLines[i]<<endl;

	if(previousEmpty && !allLines[i].empty() ){
	    vector<string> pop2=allTokens( allLines[i] ,'-');
	    if(pop2.size() != 2){
		cerr<<"ERROR: AllPairDistanceResult: Population names cannot be split "<<allLines[i]<<endl;
		exit(1);
	    }
	    
	    if(pop2[0] == "root" ){
		populationNames.push_back(pop2[1]);		
	    }
	 

	}


	previousEmpty=allLines[i].empty();

    }

    numberOfPopulations=populationNames.size();

    for(int i=0;i<numberOfPopulations;i++){
	vector<DistanceResult> * tempdist = new  vector<DistanceResult>  ();
	for(int j=0;j<numberOfPopulations;j++){
	    DistanceResult empty;
	    tempdist->push_back(empty);
	}
	distResults->push_back(tempdist);
    }

    // cout<<numberOfPopulations<<endl;
    // cout<<vectorToString(populationNames)<<endl;

    for(unsigned int k=0;k<allLines.size();k++){
	
	if(previousEmpty && !allLines[k].empty() ){
	    vector<string> pop2=allTokens( allLines[k] , '-');

	    unsigned int i = indexInVector(populationNames,pop2[0]);
	    unsigned int j = indexInVector(populationNames,pop2[1]);

	    // cout<<"i "<<i<<endl;
	    // cout<<"j "<<j<<endl;
	    // cout<<"k "<<k<<endl;

	    // cout<<distResults->size()<<endl;

	    distResults->at(i)->at(j) =  DistanceResult(allLines[k+1]+"\n"+allLines[k+2]+"\n"+allLines[k+3]+"\n"+allLines[k+4]+"\n"+allLines[k+5]+"\n"+allLines[k+6]);
	    distResults->at(j)->at(i) =  DistanceResult(allLines[k+1]+"\n"+allLines[k+2]+"\n"+allLines[k+3]+"\n"+allLines[k+4]+"\n"+allLines[k+5]+"\n"+allLines[k+6]);
	    //cout<<distResults->at(i)->at(j)<<endl;
	}

	previousEmpty=allLines[k].empty();

    }

}


int AllPairDistanceResult::getNumberOfPopulations() const{
    return numberOfPopulations;
}

vector<string> AllPairDistanceResult::getPopulationNames() const{
    return populationNames;
}

vector<string> AllPairDistanceResult::getPopulationNamesNoAnc() const{
    vector<string> namesToReturn;
    for(int i=0;i<numberOfPopulations;i++){
        if(i==1) //no anc
    	   continue;
        namesToReturn.push_back(populationNames[i]);
    }
    return namesToReturn;
}



DistanceResult AllPairDistanceResult::getDist(int pop1,int pop2) const{
    if(pop1<0 || pop1 >= numberOfPopulations){
	cerr<<"Wrong index of population in first argument: "<<pop1<<endl;
    }
    
    if(pop2<0 || pop2 >= numberOfPopulations){
	cerr<<"Wrong index of population in second argument: "<<pop2<<endl;
    }
    
    return distResults->at(pop1)->at(pop2);
}



//returns a numberOfPopulations-1 by numberOfPopulations-1 matrix
//
vector< vector<double> > AllPairDistanceResult::computeDNADist(string model) const {
    vector< vector<double> > distanceMatrix= vector< vector<double> > (numberOfPopulations-1, vector<double>(numberOfPopulations-1,0.0));



    if(model == "none"){

	for(int i=0;i<numberOfPopulations;i++){
	    int i_=i;
	    if(i==1) //no anc
	     	continue;
	    if(i>1)
		i_--;

	    vector<double> tempvec;

	    for(int j=0;j<numberOfPopulations;j++){
		int j_=j;
		if(j==1) //no anc
		    continue;
		if(j>1)
		    j_--;

		
		if(i<j)
		    distanceMatrix[i_][j_]=double(distResults->at(i)->at(j).all.getMutations());
		else
		    distanceMatrix[i_][j_]=double(distResults->at(j)->at(i).all.getMutations());
	    }
	}
	
    }else if(model == "transversions"){

	for(int i=0;i<numberOfPopulations;i++){
	    int i_=i;
	    if(i==1) //no anc
	     	continue;
	    if(i>1)
		i_--;

	    vector<double> tempvec;

	    for(int j=0;j<numberOfPopulations;j++){
		int j_=j;
		if(j==1) //no anc
		    continue;
		if(j>1)
		    j_--;

		
		if(i<j)
		    distanceMatrix[i_][j_]=double(distResults->at(i)->at(j).transversions.getMutations());
		else
		    distanceMatrix[i_][j_]=double(distResults->at(j)->at(i).transversions.getMutations());
	    }
	}
	
    }else if(model == "JC69"){


	for(int i=0;i<numberOfPopulations;i++){
	    int i_=i;
	    if(i==1) //no anc
		continue;
	    if(i>1)
		i_--;

	    vector<double> tempvec;

	    for(int j=0;j<numberOfPopulations;j++){
		int j_=j;
		if(j==1) //no anc
		    continue;
		if(j>1)
		    j_--;

		if(i==j){
		    distanceMatrix[i_][j_]=0.0;
		    continue;
		}

		if(i<j){	  

		    double P = double(distResults->at(i)->at(j).all.getMutations())   / double(distResults->at(i)->at(j).all.getIndent() );
		    // if(P>1){
		    // 	cerr<<"Cannot compute Jukes Cantor metric when there are more mutations than identical bases"<<endl;
		    // 	exit(1);
		    // }
		    double inLog=(1.0-P*double(4.0)/double(3.0));


		    if (inLog != 0.0 ){
			distanceMatrix[i_][j_] = -0.75*log(inLog);
		    }else{
			distanceMatrix[i_][j_] = 0.0;
		    }
		}else{
		    double P = double(distResults->at(j)->at(i).all.getMutations())   / double( distResults->at(j)->at(i).all.getIndent() );
		    double inLog=(1.0-P*double(4.0)/double(3.0));

		    if (inLog != 0.0){
			distanceMatrix[i_][j_] = -0.75*log(inLog);
		    }else{
			distanceMatrix[i_][j_] = 0.0;
		    }
		}
	    }//for each pop j
	}//for each pop i
	

    }else if(model == "K80"){
	cerr<<"to check"<<endl;
	exit(1);
	for(int i=0;i<numberOfPopulations;i++){
	    int i_=i;
	    if(i==1) //no anc
		continue;
	    if(i>1)
		i_--;

	    vector<double> tempvec;

	    for(int j=0;j<numberOfPopulations;j++){
		int j_=j;
		if(j==1) //no anc
		    continue;
		if(j>1)
		    j_--;

		if(i==j){
		    distanceMatrix[i_][j_]=0.0;
		    continue;
		}

		if(i<j){
		    // distanceMatrix[i_][j_]=double(distResults->at(i)->at(j).all.getMutations());

		    double P = double(distResults->at(i)->at(j).transitions.getMutations())   / double(distResults->at(i)->at(j).all.getIndent() + distResults->at(i)->at(j).all.getMutations());
		    double Q = double(distResults->at(i)->at(j).transversions.getMutations()) / double(distResults->at(i)->at(j).all.getIndent() + distResults->at(i)->at(j).all.getMutations());
		    double inLog1=(1.0-2.0*P-Q);
		    double inLog2=(1.0-2.0*Q);
		    if (inLog1 != 0.0 && inLog2 != 0.0){
			distanceMatrix[i_][j_] = -0.5*log(inLog1) -0.25*log(inLog2);
		    }else{
			distanceMatrix[i_][j_] = 0.0;
		    }
		}else{
		    //distanceMatrix[i_][j_]=double(distResults->at(j)->at(i).all.getMutations());
		    // distanceMatrix[i_][j_]=double(distResults->at(j)->at(i).all.getMutations());

		    double P = double(distResults->at(j)->at(i).transitions.getMutations())   / double(distResults->at(j)->at(i).all.getIndent() + distResults->at(j)->at(i).all.getMutations());
		    double Q = double(distResults->at(j)->at(i).transversions.getMutations()) / double(distResults->at(j)->at(i).all.getIndent() + distResults->at(j)->at(i).all.getMutations());
		    double inLog=((-2.0*P-Q+1.0) * pow ((-2.0*Q+1.0), 0.5));
		    if (inLog != 0.0){
			distanceMatrix[i_][j_] = -0.5*log( inLog);
		    }else{
			distanceMatrix[i_][j_] = 0.0;
		    }
		}
	    }//for each pop j
	}//for each pop i


    }else if(model == "HKY85"){
	cerr<<"to implement"<<endl;
	exit(1);
	

    }else{
	cerr<<"Undefined model"<<endl;
	exit(1);
    }
	    
    return distanceMatrix;
}


unsigned int AllPairDistanceResult::indexInVector(vector<string> v,string t){
    for(unsigned int i = 0;i<v.size();i++){
	if(v[i] == t){
	    return i;
	}
    }

    cerr<<"AllPairDistanceResult::indexInVector  cannot find "<<t<<" in  "<<vectorToString(v)<<endl;
    exit(1);
}
