/*
 * mistar2vcf
 * Date: Apr-03-2016
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "MistarParser.h"

using namespace std;

int main (int argc, char *argv[]) {
    string usage=string("\t"+string(argv[0])+"  [mistar file] \n"
			"This program produces a matrix where each record for population becomes a single allele {A,C,G,T}.");

     
    if(argc != 2 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }



    MistarParser mp (argv[1]);
    AlleleRecords * record;
    cout<<"##fileformat=VCFv4.0"<<endl;
    cout<<"##Tassel=<ID=GenotypeTable,Version=5,Description=\"Reference allele is not known. The major allele was used as reference allele\">"<<endl;
    cout<<"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"<<endl;
    cout<<"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">"<<endl;
    cout<<"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">"<<endl;
    cout<<"##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">"<<endl;
    cout<<"##FORMAT=<ID=PL,Number=3,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">"<<endl;
    cout<<"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"<<endl;
    cout<<"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"<<endl;
    cout<<"##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">"<<endl;
    

    //    cout<<mp.getDefline()<<endl;
    cout<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
    //#chr\tcoord\tREF,ALT\t";
    vector<string> toprintpop;
    for(unsigned int i=0;i<mp.getPopulationsNames()->size();i++){
	 toprintpop.push_back(mp.getPopulationsNames()->at(i));
     } 
     cout<<vectorToString(toprintpop,"\t")<<endl;

     while(mp.hasData()){

	record = mp.getData();
	cout<<record->chr<<"\t";
	cout<<stringify(record->coordinate)<<"\t";
	cout<<"SI_"<<stringify(record->coordinate)<<"\t";
	cout<<record->ref<<"\t";
	cout<<record->alt<<"\t";
	cout<<"."<<"\t";
	cout<<"PASS"<<"\t";
	cout<<"."<<"\t";
	cout<<"GT"<<"\t";

	vector<string> toprint;
	for(unsigned int i=0;i<record->vectorAlleles->size();i++){
	    toprint.push_back(record->vectorAlleles->at(i).generateVCFAllele(record->ref,record->alt));
	} 
	cout<<vectorToString(toprint,"\t")<<endl;
	//	}
	//
    }

     cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}

