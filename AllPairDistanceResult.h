/*
 * AllPairDistanceResult
 * Date: May-01-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef AllPairDistanceResult_h
#define AllPairDistanceResult_h

#include <float.h>
#include <cmath>
#include <vector>
#include <fstream>
#include "DistanceResult.h"

#include "utils.h"

using namespace std;

class AllPairDistanceResult{
 private:
    //DistanceResult  distResults[numberOfPopulations][numberOfPopulations];
    int numberOfPopulations;
    vector<string> populationNames;

 public:
    vector< vector<DistanceResult> * > *  distResults;

    AllPairDistanceResult(int numberTaxon,vector<string> _populationNames);
    AllPairDistanceResult(const AllPairDistanceResult & other);
    AllPairDistanceResult(const string stringForm);

    ~AllPairDistanceResult();
    AllPairDistanceResult & operator= (const AllPairDistanceResult & other);
    int getNumberOfPopulations() const;
    vector< vector<double> > computeDNADist(string model="none") const;
    vector<string> getPopulationNames() const;
    vector<string> getPopulationNamesNoAnc() const;

    DistanceResult getDist(int pop1,int pop2) const;

    AllPairDistanceResult &  operator+=(const AllPairDistanceResult & other){

	if(this->numberOfPopulations != other.numberOfPopulations){
	    cerr<<"Cannot add two AllPairDistanceResult if the number of populations differ "<<this->numberOfPopulations<<" vs "<<other.numberOfPopulations<<endl;
	    exit(1);
	}

	if(this->populationNames != other.populationNames){
	    cerr<<"Cannot add two AllPairDistanceResult if the number of populations differ "<<this->numberOfPopulations<<" vs "<<other.numberOfPopulations<<endl;
	    exit(1);
	}

	for(int i=0;i<numberOfPopulations;i++){
	    for(int j=0;j<numberOfPopulations;j++){
		/* cout<<"before "<<this->distResults->at(i)->at(j)<<endl; */
		/* cout<<other.distResults->at(i)->at(j)<<endl; */
		this->distResults->at(i)->at(j)+=other.distResults->at(i)->at(j);
		/* cout<<"after  "<<this->distResults->at(i)->at(j)<<endl; */
		//exit(1);
	    }
	}

	return *this;
    }

    unsigned int indexInVector(vector<string> v,string t);

    friend ostream& operator<<(ostream& os, const AllPairDistanceResult & adr){
	/* os<<adr.numberOfPopulations<<endl; */
	for(int i=0;i<adr.numberOfPopulations;i++){
	    if(i==1) //no anc
		continue;

	    for(int j=0;j<adr.numberOfPopulations;j++){	       
		if(j==1) //no anc
		    continue;

		if(j>=(i+1)){
		    os<<adr.populationNames[i]<<"-"<<adr.populationNames[j]<<endl;
		    os<<adr.distResults->at(i)->at(j)<<endl;
		}
	    }
	}
       
	return os;
    }


};
#endif
