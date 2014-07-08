/*
 * VecAllPairDistanceResult
 * Date: Jun-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef VecAllPairDistanceResult_h
#define VecAllPairDistanceResult_h

#include "AllPairDistanceResult.h"

using namespace std;

class VecAllPairDistanceResult{
private:
    vector<AllPairDistanceResult * > vectorAPDR;
    static string delim_;
    static string delim;

public:
    VecAllPairDistanceResult();
    VecAllPairDistanceResult(const string & vapdr );

    VecAllPairDistanceResult(const VecAllPairDistanceResult & other);
    ~VecAllPairDistanceResult();
    VecAllPairDistanceResult & operator= (const VecAllPairDistanceResult & other);

    unsigned int size() const;
    void  push_back(AllPairDistanceResult * toadd);
    const AllPairDistanceResult *  at(unsigned int i) const;

    friend std::ostream & operator<<(std::ostream & os, const VecAllPairDistanceResult & vcapdr);

};


#endif
