==========================================================
  MISTARtools: a suite of utilities for the management of allele frequency information
==========================================================

QUESTIONS :
   gabriel [dot] reno [ at sign ] gmail.com


About
----------------------

MISTARtools is a suite of utilities to convert various file formats (VCF,BAM,Affymetrix) into allele frequency matrices. MISTARtools also provides utilities to filter, combine and compute summary statistics on those matrices. These matrices can be exported to various formats for population genetics applications (treemix,fasta,EIGENSTRAT,PLINK).

Downloading:
----------------------

Go to https://github.com/grenaud/mistartools and either:

1) Download ZIP 

or

2) Do a "git clone --recursive https://github.com/grenaud/mistartools.git"

Installation
----------------------

1) make sure you have "cmake" and "git" installed, check for it by typing " git --version" and "cmake --version"

2) Make sure you have Armadillo C++ library installed. On Ubuntu, just type:

sudo apt-get install libarmadillo2 libarmadillo-dev

3) Build the submodules and main code by typing :
    make

4) Although it is not necessary, we recommend having the "bgzip" and "tabix" commands as part of your path to facilitate processing. 

Documentation
-----------------

The documentation is found here:

     doc/reference.pdf


Example of usage
-----------------

We have 4 VCF files as testData:

testData/Altai.vcf.gz
testData/Denisova.vcf.gz
testData/French.vcf.gz
testData/Yoruba.vcf.gz

For the Altai Neandertal, the Denisova, a French individual and a Yoruba individual, all high coverage genomes. 

- Convert the VCF files to MISTARtools files:

./vcf2mistar testData/Altai.vcf.gz    AltaiNeandertal testData/chr21.epo.gz |bgzip -c > testData/Altai.mst.gz 
./vcf2mistar testData/Denisova.vcf.gz Denisova        testData/chr21.epo.gz |bgzip -c > testData/Denisova.mst.gz 
./vcf2mistar testData/French.vcf.gz   French          testData/chr21.epo.gz |bgzip -c > testData/French.mst.gz 
./vcf2mistar testData/Yoruba.vcf.gz   Yoruba          testData/chr21.epo.gz |bgzip -c > testData/Yoruba.mst.gz 


- Tabix index them:

tabix -s 1 -b 2 -e 2 testData/Altai.mst.gz 
tabix -s 1 -b 2 -e 2 testData/Denisova.mst.gz 
tabix -s 1 -b 2 -e 2 testData/French.mst.gz 
tabix -s 1 -b 2 -e 2 testData/Yoruba.mst.gz 

- Create the intersection:

./mistarintersect  testData/{Altai,Denisova,French,Yoruba}.mst.gz |bgzip -c > testData/all.mst.gz


- Merge the modern humans and archaic as one population:

./mistarmeld testData/all.mst.gz   "AltaiNeandertal,Denisova" "Archaics"  |./mistarmeld /dev/stdin   "French,Yoruba" "Modern"|bgzip -c > testData/all.merged.mst.gz


- Visualize sites where the archaics and modern differ:

./mistarfilter  znosharing testData/all.merged.mst.gz "Archaics"  "Modern"


- Visualize sites where the archaics and modern differ and the archaic is ancestral and the modern humans are derived:

./mistarfilter  znosharing testData/all.merged.mst.gz "Archaics"  "Modern"  | ./mistarfilter sharing /dev/stdin  "root" "Archaics"


- Visualize sites where the archaics and modern differ and the archaic is derived and the modern humans are ancestral:

./mistarfilter  znosharing testData/all.merged.mst.gz "Archaics"  "Modern"  | ./mistarfilter sharing /dev/stdin  "root" "Modern"


- Build a neighbor-joining tree:

./mistarcompute nj  --model none testData/all.mst.gz  > testData/all.nw 


- Export to treemix:

./mistar2treemix  testData/all.mst.gz  |gzip > testData/all.treemix.gz


