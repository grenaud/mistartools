/*
 * AvgCoaResult
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef AvgCoaResult_h
#define AvgCoaResult_h

using namespace std;

#include "AlleleCounter.h"

class AvgCoaResult{
 private:
 public:
    //This is everything
    AlleleCounter all;
    //This is without the ones marked as CpG
    AlleleCounter noCpg;
    //This is only with the ones marked as CpG
    AlleleCounter onlyCpg;
    //This excludes the following cases:
    // S = C, R or A = T
    // S = T, R or A = C
    // S = A, R or A = G
    // S = G, R or A = A    
    AlleleCounter transversions;


    AlleleCounter transitions;


    //This excludes the following cases:
    // S = T, R or A = C
    // S = A, R or A = G
    AlleleCounter noDamage;

    AvgCoaResult();
    ~AvgCoaResult();
    string getHeader();

    friend ostream& operator<<(ostream& os, const AvgCoaResult & ct){
	os<<ct.all<<"\t"
	  <<ct.onlyCpg<<"\t"
	  <<ct.noCpg<<"\t"
	  <<ct.transitions<<"\t"
	  <<ct.transversions<<"\t"
	  <<ct.noDamage;
	return os;
    }
};
#endif
