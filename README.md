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

sudo apt-get install libarmadillo2
sudo apt-get install libarmadillo-dev

3) Build the submodules and main code by typing :
    make

Documentation
-----------------

The documentation is found here:

     doc/reference.pdf

