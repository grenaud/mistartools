/*
 * NewickParser
 * Date: Mar-20-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "NewickParser.h"



NewickParser::NewickParser(){

}

NewickParser::~NewickParser(){

}




Tree * NewickParser::parseFile(string filename){
    // cout<<filename<<endl;
    string line;
    igzstream myFile;
    myFile.open(filename.c_str(), ios::in);
    NodeTree * root=0;
    if (myFile.good()){
	while ( getline (myFile,line)){
	    // cout<<line<<endl;
	    if(line[line.size()-1] == ';')
		line=line.substr(0,line.size()-1)+":0";
	    root=parseNewick(line);
	    break;
	}
	myFile.close();
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	exit(1);
    }

    //cout<<"root "<<*root<<endl;

    Tree * newtree=new Tree(root);   
    return newtree;
}

