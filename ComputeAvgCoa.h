/*
 * ComputeAvgCoa
 * Date: Aug-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef ComputeAvgCoa_h
#define ComputeAvgCoa_h

#include "FilterVCF.h"
#include "AvgCoaResult.h"
#include "SimpleVCF.h"
#include "ReadTabix.h"
#include "VCFreader.h"
#include "BAMTABLEreader.h"
//#include "AlleleCounter.h"
#include "GenomicRange.h"
#include "ComputeAvgCoa_core.h"


using namespace std;

int computeAvgCoa(string refereVCF,string refereVCFidx,string sampleVCF,string sampleVCFidx,bool   useBedFileRegions,vector<GenomicRange> *  bedRegionsToFilter,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,string chrName, int startChrCoord,  int endChrCoord,SetVCFFilters * filtersVCFREF,SetVCFFilters * filtersVCFSMP,int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,int minPLdiffind,bool maximizeDiv);


int computeAvgCoa(string refereVCF,string refereVCFidx,string sampleVCF,string sampleVCFidx,bool   useBedFileRegions,vector<GenomicRange> *  bedRegionsToFilter,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,GenomicRange grc,SetVCFFilters * filtersVCFREF,SetVCFFilters * filtersVCFSMP,int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,int minPLdiffind,bool maximizeDiv);

#endif
