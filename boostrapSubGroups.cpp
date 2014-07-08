/*
 * testNewickParser
 * Date: Mar-20-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>
#include "NewickParser.h"
#include "NodeTree.h"
#include "Tree.h"

#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {
    

    const string usage=string(string(argv[0])+" [reference tree] [bootstrap tree 1] [bootstrap tree 2] ...");

    if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
    	cout<<"Usage:"<<endl;
    	cout<<""<<endl;
    	cout<<usage<<endl;
    	return 1;
    }

    

    //    cout<<argv[1]<<endl;
    NewickParser np;

    Tree * reference=np.parseFile(argv[1]);
    vector<Tree *> bootReplicates;
    vector< set<string> > subtrees = reference->returnAllSubgroup();
    vector<int> count (subtrees.size(),0);


    for(int i=2;i<(argc);i++){
	Tree * currentReplicate=np.parseFile(argv[i]);
	//bootReplicates->push_back(currentReplicate);
	for(unsigned int i=0;i<subtrees.size();i++){
	    if(currentReplicate->hasSubGroup( subtrees[i] )){
		count[i]++;
	    }
	}
	delete(currentReplicate);
    }

        
    for(unsigned int i=0;i<subtrees.size();i++){

	if(subtrees[i].size() > 1){
	    cout<<"{"<<iteratorToString(subtrees[i])<<"}\t"<<count[i]<<endl;
	}
    }
    
    delete(reference);
    
    return 0;
}

