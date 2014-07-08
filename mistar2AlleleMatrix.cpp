/*
 * mistar2AlleleMatrix
 * Date: Mar-28-2013 
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
    //    cout<<mp.getDefline()<<endl;
    cout<<"#chr\tcoord\tREF,ALT\t";
    vector<string> toprintpop;
    for(unsigned int i=0;i<mp.getPopulationsNames()->size();i++){
	 toprintpop.push_back(mp.getPopulationsNames()->at(i));
     } 
     cout<<vectorToString(toprintpop,"\t")<<endl;

     while(mp.hasData()){
	//cout<<"data"<<endl;
	record = mp.getData();
	

	//if(record->everyNonChimpAncRecordNonNull()){
	//     //	   
	//	    cout<<*record<<endl;
	cout<<record->chr<<"\t";
	cout<<stringify(record->coordinate)<<"\t";
	cout<<record->ref<<",";
	cout<<record->alt<<"\t";
	vector<char> toprint;
	for(unsigned int i=0;i<record->vectorAlleles->size();i++){
	    toprint.push_back(record->vectorAlleles->at(i).generateRandomAlleleBias(record->ref,record->alt));
	} 
	cout<<vectorToString(toprint,"\t")<<endl;
	//	}
	//
    }

     cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}

