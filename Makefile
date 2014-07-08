CXX      = g++ #-g -pg


LIBGAB   = libgab/
LIBTABIX = tabix/
BAMTOOLS = bamtools/

LDFLAGS  =   ReadTabix.o  tabix/libtabix.a /usr/lib/libarmadillo.so.2    -lz
CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB}/   -I${LIBTABIX}   -I${LIBGAB}/gzstream/  -I${BAMTOOLS}/include/  -I${BAMTOOLS}/src/  -c

all: tabix/libtabix.a bamtools/lib/libbamtools.a libgab/utils.o ReadTabix.o testNewickParser boostrapSubGroups mistar2treemix mistarcompute mistarcat  mistarintersect mistarmeld mistarunion testMistar mistarfilter mistar2AlleleMatrix epo2mistar mistaruniq ms2nj mistar2binary closestMistar replaceAncestor mistar2nexus ms2mistar mistar2bed mistarfreqSpec  testMistarTabix mistarcompress mistardecompress mpileup2mistar mistar2EIGENSTRAT mistar2BinaryPLINK mistar2segtor mistar2fasta mistarRenamePop usePopAsRootAnc mistarstats axt2mistar mistarSubsample bam2mistar vcfcompute testRandomCoordGenerator testComputation generateCoords testMS avgCoaMS testReadTabix testVCF mergeBAMTable filterVCF testReadFastq vcf2mistar bamtable2mistar affyVCF2mistar 23andme2mistar testMultiReadTabix vcfMulti2mistar


%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@


bamtools/src/api/BamAlignment.h:
	rm -rf bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git


bamtools/lib/libbamtools.a: bamtools/src/api/BamAlignment.h
	cd bamtools/ && mkdir build/  && cd build/ && cmake .. && make && cd ../..


libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git

libgab/utils.o: bamtools/lib/libbamtools.a  libgab/utils.h
	make -C libgab

tabix/tabix.h:
	rm -rf tabix/
	git clone --recursive https://github.com/samtools/tabix.git

tabix/libtabix.a: tabix/tabix.h
	make -C tabix




#testGMM mistar




vcfcompute: vcfcompute.o libgab/utils.o    VCFreader.o SimpleVCF.o  CoreVCF.o  BAMTABLEreader.o  BAMTableObj.o tabix/libtabix.a  AlleleCounter.o ComputeAvgCoa.o RandomGenomicCoord.o GenomicRange.o ComputeFracAnc.o ComputeFracHetero.o  FilterVCF.o  GenomicWindows.o MultiVCFParser.o Dstats.o DstatResult.o DstatCounter.o Dstat_core.o SetVCFFilters.o AvgCoaResult.o ComputeAvgCoa_core.o libgab//gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS)  $(LDFLAGS)


testMS:	testMS.o libgab/utils.o MSParser.o MSobject.o libgab//gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 

avgCoaMS:	avgCoaMS.o libgab/utils.o MSParser.o MSobject.o libgab//gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 


generateCoords:	generateCoords.o libgab/utils.o  GenomicWindows.o RandomGenomicCoord.o GenomicRange.o libgab//gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testComputation:	testComputation.o libgab/utils.o ComputeFracAnc.o ComputeFracHetero.o  VCFreader.o CoreVCF.o SimpleVCF.o BAMTABLEreader.o BAMTableObj.o libgab//gzstream/libgzstream.a tabix/libtabix.a GenomicRange.o FilterVCF.o  SetVCFFilters.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

testRandomCoordGenerator:	testRandomCoordGenerator.o libgab/utils.o   RandomGenomicCoord.o  GenomicRange.o GenomicWindows.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 


#${LIBMLPACK}/build/lib/libmlpack.so  

axt2mistar.o:	axt2mistar.cpp
	${CXX} ${CXXFLAGS} axt2mistar.cpp


axt2mistar:	axt2mistar.o libgab/utils.o  libgab//gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


testGMM: testGMM.o ComputeGMM.o mistarArmadillo.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistar: mistar.o  libgab/utils.o  MistarParser.o libgab//gzstream/gzstream.o AllelePairCounter.o AlleleCounter.o ComputeAvgCoa_core.o SingleAllele.o AlleleRecords.o ComputeGMM.o MistarPairwiseDiff.o AllPairDistanceResult.o VecAllPairDistanceResult.o DistanceResult.o mistarOperations.o mistarArmadillo.o NjTree.o Tree.o NodeTree.o UnrootedNode.o GenomicRange.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistarcompute:	mistarcompute.o libgab/utils.o MistarPairwiseAvgCoa.o  MistarPairwiseDiff.o MistarParser.o libgab//gzstream/gzstream.o AllelePairCounter.o DistanceResult.o AvgCoaResult.o AlleleCounter.o ComputeAvgCoa_core.o DstatResult.o DstatCounter.o MistarDstats.o Dstat_core.o SingleAllele.o AlleleRecords.o NjTree.o Tree.o NodeTree.o UnrootedNode.o AllPairDistanceResult.o  mistarArmadillo.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

ms2nj:	ms2nj.o	 libgab/utils.o libgab//gzstream/gzstream.o MSParser.o MSobject.o NjTree.o Tree.o NodeTree.o UnrootedNode.o AlleleCounter.o AllelePairCounter.o  mistarArmadillo.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

ms2mistar:	ms2mistar.o	 libgab/utils.o libgab//gzstream/gzstream.o MSParser.o MSobject.o NjTree.o Tree.o NodeTree.o UnrootedNode.o AlleleCounter.o AllelePairCounter.o SingleAllele.o mistarArmadillo.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mpileup2mistar:	mpileup2mistar.o	 libgab/utils.o libgab//gzstream/gzstream.o MSParser.o MSobject.o NjTree.o Tree.o NodeTree.o UnrootedNode.o AlleleCounter.o AllelePairCounter.o SingleAllele.o mistarArmadillo.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistar2treemix:	mistar2treemix.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

closestMistar:	closestMistar.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 


bam2mistar:	bam2mistar.o libgab/utils.o  libgab//gzstream/gzstream.o   bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o bamtools/lib/libbamtools.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarmeld:	mistarmeld.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarintersect:	mistarintersect.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarfilter:	mistarfilter.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarfreqSpec:	mistarfreqSpec.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2bed:	mistar2bed.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2binary:	mistar2binary.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2segtor:	mistar2segtor.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarSubsample:	mistarSubsample.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2EIGENSTRAT:	mistar2EIGENSTRAT.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2BinaryPLINK:	mistar2BinaryPLINK.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarcompress:	mistarcompress.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistardecompress:	mistardecompress.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistar2nexus:	mistar2nexus.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarunion:	mistarunion.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

replaceAncestor:	replaceAncestor.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistaruniq:	mistaruniq.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o mistarOperations.o GenomicRange.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

mistarcat:	mistarcat.o libgab/utils.o   MistarParser.o libgab//gzstream/libgzstream.a SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testNewickParser:	testNewickParser.o libgab/utils.o Tree.o  NodeTree.o UnrootedNode.o NewickParser.o  libgab//gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testMistar:	testMistar.o libgab/utils.o   MistarParser.o libgab//gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testMistarTabix:	testMistarTabix.o libgab/utils.o   MistarParser.o libgab//gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistar2AlleleMatrix:	mistar2AlleleMatrix.o libgab/utils.o   MistarParser.o libgab//gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

epo2mistar:	epo2mistar.o libgab/utils.o   MistarParser.o libgab//gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistar2AlleleMatrix:	mistar2AlleleMatrix.o   MistarParser.o libgab/utils.o  libgab//gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistar2fasta:	mistar2fasta.o libgab/utils.o  MistarParser.o libgab//gzstream/libgzstream.a  SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

boostrapSubGroups:	boostrapSubGroups.o libgab/utils.o Tree.o  NodeTree.o UnrootedNode.o NewickParser.o libgab//gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistarRenamePop: mistarRenamePop.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

usePopAsRootAnc: usePopAsRootAnc.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

mistarstats: mistarstats.o libgab/utils.o MistarParser.o libgab//gzstream/gzstream.o SingleAllele.o AlleleRecords.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testVCF:	testVCF.o libgab/utils.o SimpleVCF.o CoreVCF.o FilterVCF.o BAMTableObj.o SetVCFFilters.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testReadTabix:	testReadTabix.o libgab/utils.o  tabix/libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o  libgab//gzstream/gzstream.o 
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testMultiReadTabix: testMultiReadTabix.o MultiVCFreader.o libgab/utils.o  tabix/libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o  libgab//gzstream/gzstream.o 
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

mergeBAMTable:	mergeBAMTable.o libgab/utils.o BAMTableObj.o BAMTABLEreader.o    tabix/libtabix.a libgab//gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

filterVCF:	filterVCF.o libgab/utils.o    tabix/libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o FilterVCF.o SetVCFFilters.o libgab//gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

vcf2mistar:	vcf2mistar.o libgab/utils.o   tabix/libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o FilterVCF.o SetVCFFilters.o libgab//gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

vcfMulti2mistar: 	vcfMulti2mistar.o libgab/utils.o  MultiVCFreader.o   tabix/libtabix.a SimpleVCF.o CoreVCF.o   BAMTableObj.o BAMTABLEreader.o FilterVCF.o SetVCFFilters.o libgab//gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)


affyVCF2mistar:	affyVCF2mistar.o libgab/utils.o   tabix/libtabix.a  libgab//gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

23andme2mistar:	23andme2mistar.o libgab/utils.o   tabix/libtabix.a  libgab//gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

bamtable2mistar:	bamtable2mistar.o libgab/utils.o   tabix/libtabix.a   BAMTableObj.o BAMTABLEreader.o  libgab/gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)


clean :
	rm -f  *.o  testNewickParser boostrapSubGroups mistar2treemix mistarcompute mistarcat  mistarintersect mistarmeld mistarunion testMistar mistarfilter mistar2AlleleMatrix epo2mistar mistaruniq ms2nj mistar2binary closestMistar replaceAncestor mistar2nexus ms2mistar mistar2bed mistarfreqSpec  testMistarTabix mistarcompress mistardecompress mpileup2mistar mistar2EIGENSTRAT mistar2BinaryPLINK mistar2segtor mistar2fasta mistarRenamePop usePopAsRootAnc mistarstats axt2mistar mistarSubsample bam2mistar vcfcompute testRandomCoordGenerator testComputation generateCoords testMS avgCoaMS testReadTabix testVCF mergeBAMTable filterVCF testReadFastq vcf2mistar bamtable2mistar affyVCF2mistar 23andme2mistar testMultiReadTabix vcfMulti2mistar
