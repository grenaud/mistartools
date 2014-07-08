/*
 * mistarcompute
 * Date: Jan-26-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "MistarPairwiseDiff.h"
#include "MistarPairwiseAvgCoa.h"
#include "MistarDstats.h"

using namespace std;

int main (int argc, char *argv[]) {
    char mode; //d = divergence, p=pairdiff, s=d-stats, t= nj tree
    string filenameInput;
    string dnaDistMode;
    bool allowUndefined=false;

    string usage=string(""+string(argv[0])+"  [mode]"+
			"\nThis program has the following modes:\n\n"+
			"\t\t"+"paircoacompute  To compute pairwise average coalescence\n"+
			"\t\t"+"pairdiff        To compute pairwise nucleotide differences\n"+
			"\t\t"+"nj              Compute neighbor-joining tree\n"+
			"\t\t"+"dstat           To compute triple-wise D-statistics\n");

    string usageModel=string("\t\t")+"--model [model]"+"\t"+"Use this model for DNA distance\n"+"\n"+
	"\t\t\t"+"none"+"\t\t"+"all mutations with equal footing\n"+
	"\t\t\t"+"transversions"+"\t"+"Just transversions\n"+
	"\t\t\t"+"JC69"+"\t\t"+"Jukes Cantor 1969\n"+
	"\t\t\t"+"K80"+"\t\t"+"Kimura 1980\n"+
	"\t\t\t"+"HKY85"+"\t\t"+"Hasegawa, Kishino and Yano 1985\n";

    string optionAvgCoa=string("\t\t")+"--undef"+"\t"+"Ignore undefined sites (cannot compare between pairs)\n";

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    if( string(argv[1]) == "paircoacompute"){
	mode='d';	
	usage=string(string(argv[0])+" "+string(argv[1])+" [options] <mistar file>\n\n");
    }else{
	if( string(argv[1]) == "pairdiff"){
	    mode='p';
	    usage=string(string(argv[0])+" "+string(argv[1])+" [options] <mistar file>\n\n");
	}else{
	    if( string(argv[1]) == "dstat"){
		mode='s';
		usage=string(string(argv[0])+" "+string(argv[1])+" [options] <mistar file>\n\n");
	    }else{
		if( string(argv[1]) == "nj"){
		    mode='t';
		    usage=string(string(argv[0])+" "+string(argv[1])+" [options] <mistar file>\n\n"+usageModel);
		}else{
		    cerr << "Invalid mode\nUsage:\n"<<usage<<endl;
		    return 1; 
		}
	    }
	}
    }


    if(argc == 2 ||
       (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    //starts at two for the sub-program name
    for(int i=2;i<(argc-1);i++){ 
	//cout<<string(argv[i])<<endl;
	if(mode == 't'){
	    if(string(argv[i]) == "--model"){
		dnaDistMode=string(argv[i+1]);
		i++;
		continue;
	    }
	}

	if(mode == 'd'){
	    if(string(argv[i]) == "--undef"){

		allowUndefined=true;
		continue;
	    }
	}

	cerr<<"Unknown command line argument  "<<argv[i]<<endl;
	return 1;
    }

    filenameInput=string(argv[argc-1]);

    
    //BEGIN summary of parameters used
    cerr<<"Parameters used:"<<endl;
    cerr<<"Date: "<<getDateString()<<"\t"<<getTimeString()<<endl;
    cerr<<"Github version: "<<returnGitHubVersion(argv[0],"")<<endl;
    AllPairDistanceResult * apdr;
    vector< vector<double> > dist;

    switch(mode){

    case 'd': //divergence
	
	pairwiseAvgCoa(filenameInput,allowUndefined);
	break;

    case 'p': //differences

	cout<<*pairwiseDifferences(filenameInput)<<endl;
	break;

    case 't': //tree

	//pairwiseDifferences(filenameInput,true);
	apdr=pairwiseDifferences(filenameInput);
	// cout<<*apdr<<endl;
	dist=apdr->computeDNADist(dnaDistMode);
	// for(unsigned i=0;i<dist.size();i++){
	//     cout<<vectorToString(dist[i],"\t")<<endl;
	// }

	cout<<*neighborJoinFromDist(dist,
				    apdr->getPopulationNamesNoAnc(),
				    apdr->getNumberOfPopulations()-1)<<endl;
	break;

    case 's': //d-stats

	dstatsMistar(filenameInput);
	break;

    default:
	cerr << "Unknown mode "<<mode<<endl;
	return 1;       
	break;	 
    }

	

	
    cerr<<"Program terminated gracefully"<<endl;

    

    return 0;
}

