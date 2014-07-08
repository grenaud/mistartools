/*
 * VecAllPairDistanceResult
 * Date: Jun-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "VecAllPairDistanceResult.h"

string VecAllPairDistanceResult::delim_="/--------------------------------------------------/";
string VecAllPairDistanceResult::delim="\n"+delim_+"\n";

VecAllPairDistanceResult::VecAllPairDistanceResult(){
    // delim="/---------------------------/";
}

VecAllPairDistanceResult::VecAllPairDistanceResult(const string & vapdr){
    //cerr<<vapdr<<endl;
    
    istringstream f (vapdr);
    string line; 
    string accumulatorLine="";
    while (getline(f, line)) {
	//cout<<"line "<<line<<endl;
	if(line == delim_){//found new record
	    AllPairDistanceResult * toadd = new  AllPairDistanceResult(accumulatorLine);
	    push_back(toadd);
	    //cout<<"toadd "<<*toadd<<endl;
	    accumulatorLine="";	    
	}else{
	    accumulatorLine+=line+"\n";
	}

    }

    if(!accumulatorLine.empty()){
	AllPairDistanceResult * toadd = new  AllPairDistanceResult(accumulatorLine);
	push_back(toadd);
    }
}



VecAllPairDistanceResult::~VecAllPairDistanceResult(){
}


unsigned int VecAllPairDistanceResult::size() const{
    return vectorAPDR.size();
}


void  VecAllPairDistanceResult::push_back(AllPairDistanceResult * toadd){
    vectorAPDR.push_back(toadd);
}


const AllPairDistanceResult *  VecAllPairDistanceResult::at(unsigned int i) const{
    return vectorAPDR[i];
}

std::ostream & operator<<(std::ostream & os, const VecAllPairDistanceResult & vcapdr){
    if(vcapdr.vectorAPDR.size() == 0){
	os<<"";
    }

    for(unsigned int i=0;i<((vcapdr.vectorAPDR.size())-1);i++){
	os<<*(vcapdr.vectorAPDR[i])<<vcapdr.delim;
    }

    os<<*(vcapdr.vectorAPDR[ vcapdr.vectorAPDR.size() -1 ]);

    return os;

}
