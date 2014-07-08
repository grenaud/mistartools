CXX      = g++ #-g -pg
LIBGAB   = /home/gabriel_renaud/lib/
LIBTABIX = /home/gabriel_renaud/Software/tabix-0.2.6/
BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2/



LDFLAGS  =   ReadTabix.o  ${LIBTABIX}libtabix.a /usr/lib/libarmadillo.so.2    -lz

all: ReadTabix.o testNewickParser boostrapSubGroups mistar2treemix mistarcompute mistarcat  mistarintersect mistarmeld mistarunion testMistar mistarfilter mistar2AlleleMatrix epo2mistar mistaruniq ms2nj mistar2binary closestMistar replaceAncestor mistar2nexus ms2mistar mistar2bed mistarfreqSpec  testMistarTabix mistarcompress mistardecompress mpileup2mistar mistar2EIGENSTRAT mistar2BinaryPLINK mistar2segtor mistar2fasta mistarRenamePop usePopAsRootAnc mistarstats axt2mistar mistarSubsample bam2mistar vcfcompute testRandomCoordGenerator testComputation generateCoords testMS avgCoaMS testReadTabix testVCF mergeBAMTable filterVCF testReadFastq vcf2mistar bamtable2mistar affyVCF2mistar 23andme2mistar testMultiReadTabix vcfMulti2mistar

#testGMM mistar


%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@


vcfcompute: vcfcompute.o ${LIBGAB}utils.o    VCFreader.o SimpleVCF.o  CoreVCF.o  BAMTABLEreader.o  BAMTableObj.o /home/gabriel_renaud/Software/tabix-0.2.6/libtabix.a  AlleleCounter.o ComputeAvgCoa.o RandomGenomicCoord.o GenomicRange.o ComputeFracAnc.o ComputeFracHetero.o  FilterVCF.o  GenomicWindows.o MultiVCFParser.o Dstats.o DstatResult.o DstatCounter.o Dstat_core.o SetVCFFilters.o AvgCoaResult.o ComputeAvgCoa_core.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS)  $(LDFLAGS)


testMS:	testMS.o ${LIBGAB}utils.o MSParser.o MSobject.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 

avgCoaMS:	avgCoaMS.o ${LIBGAB}utils.o MSParser.o MSobject.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 


generateCoords:	generateCoords.o ${LIBGAB}utils.o  GenomicWindows.o RandomGenomicCoord.o GenomicRange.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testComputation:	testComputation.o ${LIBGAB}utils.o ComputeFracAnc.o ComputeFracHetero.o  VCFreader.o CoreVCF.o SimpleVCF.o BAMTABLEreader.o BAMTableObj.o ${LIBGAB}/gzstream/libgzstream.a /home/gabriel_renaud/Software/tabix-0.2.6/libtabix.a GenomicRange.o FilterVCF.o  SetVCFFilters.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

testRandomCoordGenerator:	testRandomCoordGenerator.o ${LIBGAB}utils.o   RandomGenomicCoord.o  GenomicRange.o GenomicWindows.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 


#${LIBMLPACK}/build/lib/libmlpack.so  

axt2mistar.o:	axt2mistar.cpp
	${CXX} ${CXXFLAGS} axt2mistar.cpp


axt2mistar:	axt2mistar.o ${LIBGAB}utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


testGMM: testGMM.o ComputeGMM.o mistarArmadillo.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistar: mistar.o  ${LIBGAB}utils.o  MistarParser.o ${LIBGAB}/gzstream/gzstream.o AllelePairCounter.o AlleleCounter.o ComputeAvgCoa_core.o SingleAllele.o AlleleRecords.o ComputeGMM.o MistarPairwiseDiff.o AllPairDistanceResult.o VecAllPairDistanceResult.o DistanceResult.o mistarOperations.o mistarArmadillo.o NjTree.o Tree.o NodeTree.o UnrootedNode.o GenomicRange.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistarcompute:	mistarcompute.o ${LIBGAB}utils.o MistarPairwiseAvgCoa.o  MistarPairwiseDiff.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o AllelePairCounter.o DistanceResult.o AvgCoaResult.o AlleleCounter.o ComputeAvgCoa_core.o DstatResult.o DstatCounter.o MistarDstats.o Dstat_core.o SingleAllele.o AlleleRecords.o NjTree.o Tree.o NodeTree.o UnrootedNode.o AllPairDistanceResult.o  mistarArmadillo.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

ms2nj:	ms2nj.o	 ${LIBGAB}utils.o ${LIBGAB}/gzstream/gzstream.o MSParser.o MSobject.o NjTree.o Tree.o NodeTree.o UnrootedNode.o AlleleCounter.o AllelePairCounter.o  mistarArmadillo.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

ms2mistar:	ms2mistar.o	 ${LIBGAB}utils.o ${LIBGAB}/gzstream/gzstream.o MSParser.o MSobject.o NjTree.o Tree.o NodeTree.o UnrootedNode.o AlleleCounter.o AllelePairCounter.o SingleAllele.o mistarArmadillo.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mpileup2mistar:	mpileup2mistar.o	 ${LIBGAB}utils.o ${LIBGAB}/gzstream/gzstream.o MSParser.o MSobject.o NjTree.o Tree.o NodeTree.o UnrootedNode.o AlleleCounter.o AllelePairCounter.o SingleAllele.o mistarArmadillo.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistar2treemix:	mistar2treemix.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

closestMistar:	closestMistar.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 


bam2mistar:	bam2mistar.o ${LIBGAB}utils.o  ${LIBGAB}/gzstream/gzstream.o   ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o ${BAMTOOLS}/lib/libbamtools.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarmeld:	mistarmeld.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarintersect:	mistarintersect.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarfilter:	mistarfilter.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarfreqSpec:	mistarfreqSpec.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2bed:	mistar2bed.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2binary:	mistar2binary.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2segtor:	mistar2segtor.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarSubsample:	mistarSubsample.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2EIGENSTRAT:	mistar2EIGENSTRAT.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2BinaryPLINK:	mistar2BinaryPLINK.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarcompress:	mistarcompress.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistardecompress:	mistardecompress.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2nexus:	mistar2nexus.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarunion:	mistarunion.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

replaceAncestor:	replaceAncestor.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistaruniq:	mistaruniq.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarcat:	mistarcat.o ${LIBGAB}utils.o   MistarParser.o ${LIBGAB}/gzstream/libgzstream.a SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testNewickParser:	testNewickParser.o ${LIBGAB}utils.o Tree.o  NodeTree.o UnrootedNode.o NewickParser.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testMistar:	testMistar.o ${LIBGAB}utils.o   MistarParser.o ${LIBGAB}/gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testMistarTabix:	testMistarTabix.o ${LIBGAB}utils.o   MistarParser.o ${LIBGAB}/gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistar2AlleleMatrix:	mistar2AlleleMatrix.o ${LIBGAB}utils.o   MistarParser.o ${LIBGAB}/gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

epo2mistar:	epo2mistar.o ${LIBGAB}utils.o   MistarParser.o ${LIBGAB}/gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistar2AlleleMatrix:	mistar2AlleleMatrix.o   MistarParser.o ${LIBGAB}utils.o  ${LIBGAB}/gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistar2fasta:	mistar2fasta.o ${LIBGAB}utils.o  MistarParser.o ${LIBGAB}/gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

boostrapSubGroups:	boostrapSubGroups.o ${LIBGAB}utils.o Tree.o  NodeTree.o UnrootedNode.o NewickParser.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistarRenamePop: mistarRenamePop.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

usePopAsRootAnc: usePopAsRootAnc.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistarstats: mistarstats.o ${LIBGAB}utils.o MistarParser.o ${LIBGAB}/gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testVCF:	testVCF.o ${LIBGAB}utils.o SimpleVCF.o CoreVCF.o FilterVCF.o BAMTableObj.o SetVCFFilters.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testReadTabix:	testReadTabix.o ${LIBGAB}utils.o  ${LIBTABIX}libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o  ${LIBGAB}/gzstream/gzstream.o 
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testMultiReadTabix: testMultiReadTabix.o MultiVCFreader.o ${LIBGAB}utils.o  ${LIBTABIX}libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o  ${LIBGAB}/gzstream/gzstream.o 
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

mergeBAMTable:	mergeBAMTable.o ${LIBGAB}utils.o BAMTableObj.o BAMTABLEreader.o    ${LIBTABIX}libtabix.a ${LIBGAB}/gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

filterVCF:	filterVCF.o ${LIBGAB}utils.o    ${LIBTABIX}libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o FilterVCF.o SetVCFFilters.o ${LIBGAB}/gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

vcf2mistar:	vcf2mistar.o ${LIBGAB}utils.o   ${LIBTABIX}libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o FilterVCF.o SetVCFFilters.o ${LIBGAB}/gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

vcfMulti2mistar: 	vcfMulti2mistar.o ${LIBGAB}utils.o  MultiVCFreader.o   ${LIBTABIX}libtabix.a SimpleVCF.o CoreVCF.o   BAMTableObj.o BAMTABLEreader.o FilterVCF.o SetVCFFilters.o ${LIBGAB}/gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)


affyVCF2mistar:	affyVCF2mistar.o ${LIBGAB}utils.o   ${LIBTABIX}libtabix.a  ${LIBGAB}/gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

23andme2mistar:	23andme2mistar.o ${LIBGAB}utils.o   ${LIBTABIX}libtabix.a  ${LIBGAB}/gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

bamtable2mistar:	bamtable2mistar.o ${LIBGAB}utils.o   ${LIBTABIX}libtabix.a   BAMTableObj.o BAMTABLEreader.o  ${LIBGAB}/gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)


clean :
	rm -f  *.o  testNewickParser boostrapSubGroups mistar2treemix mistarcompute mistarcat  mistarintersect mistarmeld mistarunion testMistar mistarfilter mistar2AlleleMatrix epo2mistar mistaruniq ms2nj mistar2binary closestMistar replaceAncestor mistar2nexus ms2mistar mistar2bed mistarfreqSpec  testMistarTabix mistarcompress mistardecompress mpileup2mistar mistar2EIGENSTRAT mistar2BinaryPLINK mistar2segtor mistar2fasta mistarRenamePop usePopAsRootAnc mistarstats axt2mistar mistarSubsample bam2mistar vcfcompute testRandomCoordGenerator testComputation generateCoords testMS avgCoaMS testReadTabix testVCF mergeBAMTable filterVCF testReadFastq vcf2mistar bamtable2mistar affyVCF2mistar 23andme2mistar testMultiReadTabix vcfMulti2mistar
