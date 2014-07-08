/*
 * MistarPairwiseDiff
 * Date: Jan-26-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef MistarPairwiseDiff_h
#define MistarPairwiseDiff_h

#include <string>

#include "DistanceResult.h"
#include "NjTree.h"
#include "AllPairDistanceResult.h"
#include "MistarParser.h"


using namespace std;

AllPairDistanceResult *  pairwiseDifferences(string filename);
AllPairDistanceResult *  pairwiseDifferences_calc( MistarParser * mp,string chrNameRepo="NONE",unsigned int startCoordRepo=0,unsigned int endCoordRepo=0  );


#endif
