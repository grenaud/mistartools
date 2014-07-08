/*
 * NewickParser
 * Date: Mar-20-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef NewickParser_h
#define NewickParser_h

#include <stdlib.h>
#include <string>
#include <gzstream.h>
#include <set>

#include "NodeTree.h"
#include "Tree.h"

#include "utils.h"

/* #define DEBUG */


using namespace std;
class Tree;

class NewickParser{
 private:
    /* bool doubleCheckParentheses(string line); */
 public:
    NewickParser();
    NewickParser(const NewickParser & other);
    ~NewickParser();
    NewickParser & operator= (const NewickParser & other);
    Tree * parseFile(string filename);

};


inline bool doubleCheckParentheses(string line){
    int numberOfOpen=0;
    int numberOfClose=0;

    for(unsigned int i=0;i<line.size();i++){
	if(line[i] == '('){
	    numberOfOpen++;
	}
	if(line[i] == ')'){
	    numberOfClose++;
	}
    }
    return (numberOfOpen !=  numberOfClose);	
}



inline NodeTree * parseNewick(string line){

#ifdef DEBUG
    cerr<<endl<<endl<<"parseNewick(  #"<<line<<"#  )"<<endl;
#endif

    if(doubleCheckParentheses(line)){
	cerr<<"Line does not have an equal number of parentheses: "<<line<<endl;
	exit(1);
    }

    set<string> allPopulationNames;

    double distance;
    string dist="";
    unsigned int k;
    for(k=(line.size()-1);k>0;k--){
	if(line[k] != ':'){
	    dist=line[k]+dist;
	}else
	    break;
    }
    if(dist.empty()){
	cerr<<"Line does not have a distance: "<<line<<endl;
	exit(1);
    }

#ifdef DEBUG
     cerr<<"dist "<<dist<<endl;
#endif

    distance=destringify<double>(dist);
#ifdef DEBUG
     cerr<<"distance "<<distance<<endl;
#endif


    line=line.substr(0,k);


    if(line[0] == '(' && line[line.size()-1] == ')'){//contains subtree

#ifdef DEBUG
     cerr<<"sub "<<line<<endl;
#endif
	// cerr<<"sub"<<endl;
	int indexComma=-1;

	if(line[1] == '('){//left child is subtree, the comma is when the # of '(' == # of ')'
	    int numberOfOpen=1;

	    for(unsigned int i=2;i<line.size();i++){
		if(line[i] == '('){
		    numberOfOpen++;
		}
		if(line[i] == ')'){
		    numberOfOpen--;
		}
		if(numberOfOpen==0){
		    if(line[i] == ','){
			indexComma=i;
			break;
		    }
		}
	    }
	    if(indexComma == -1 ){ cerr<<"Parsing error for : "<<line<<endl; exit(1);}
	}

	//left child is not a subtree, find the next comma
	else{
	    
	    for(unsigned int i=1;i<line.size();i++){
		if(line[i] == ','){
		    indexComma=i;
		    break;
		}		
	    }

	}
	// cerr<<"indexComma "<<indexComma<<endl;

	NodeTree * nt  = new NodeTree(parseNewick( line.substr(1,indexComma-1) ),
				      parseNewick( line.substr(indexComma+1, line.size() - indexComma-2) ) ,distance);
	return nt;

    }else{
#ifdef DEBUG
	cerr<<"no sub "<< line <<endl;
#endif
	NodeTree * nt  = new NodeTree(line,distance);
	if(allPopulationNames.find(line) == allPopulationNames.end()){
	    allPopulationNames.insert(line);
	}else{
	    cerr<<"ERROR: Duplicated names for tree "<<line<<endl;
	    exit(1);
	}
	
	return nt;
    }

    
}


#endif
