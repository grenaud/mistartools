/*
 * Tree
 * Date: Mar-21-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef Tree_h
#define Tree_h

#include <sstream>  
#include <memory>
#include <set>
#include <vector>

#include "NewickParser.h"
#include "NodeTree.h"
#include "UnrootedNode.h"

using namespace std;
/* NodeTree * parseNewick(string line); */

class Tree{
private:
    NodeTree * root;
    vector< set<string> > * allSubGroups;
    bool hasLeafWithNameIntern(NodeTree * node,string nameToLookFor);

public:
    Tree(NodeTree * _root);
    Tree(const UnrootedNode * _root,const UnrootedNode * _rest,double distance);

    Tree(const Tree & other);
    ~Tree();
    Tree & operator= (const Tree & other);

    bool hasSubTree(string tree);
    bool hasSubTreeWithName(NodeTree * node,string subTree);
    bool hasLeafWithName(string nameToLookFor);
    vector< set<string> > returnAllSubgroup();
    void returnAllSubgroupNode(const NodeTree * node,vector< set<string> > * accumulateSets);
    bool hasSubGroup(set<string> tocheck);
    
    const NodeTree * getRoot() const;
    
    void getDistanceToLeave(const NodeTree * node,
			    vector<double> * distForNodes,
			    vector<string> * distNames,
			    double currentDist) const;

    void findInternalBranch(const NodeTree * node) const;
    void getDistanceUpperTree(const NodeTree * node,
			      const NodeTree * callingNode,
			      vector<double> * distForNodes,
			      vector<string> * distNames,
			      double currentDist) const;
    void addNewNodeUsingDist(string newName,vector<double> distForNode,vector<string> distNames);

    string nameRoot(){
	stringstream toreturn;
	toreturn<<root->getPopulation();
	return toreturn.str();
    }

    /*********************/
    /*   print newick   */
    /*********************/
    
    string printTreeNewick(const bool nodist=false) const{
	return printSubtreeNewick(root,nodist);
    }

    string printSubtreeNewick(const NodeTree * node,const bool nodist=false) const{

	stringstream toreturn;
	if(node->isLeaf()){
	    if(nodist){
		toreturn<<node->getPopulation();
	    }else{
		toreturn<<*node;       
	    }
	}else{

	    toreturn<<"("<<printSubtreeNewick(node->getLeftChild(),nodist)<<","<<printSubtreeNewick(node->getRightChild(),nodist)<<")";
	    if(!nodist){
		toreturn<<":"<<node->getDistance();
	    }
	}
	/* toreturn<<";"; */
	return toreturn.str();

    }


    /*********************************/
    /*   print console "friendly"   */
    /********************************/

    string printTree(const bool nodist=false) const{
	return printSubtree(root,0,nodist);
    }
    string printSubtree(const NodeTree * node,const int depth=0,const bool nodist=false) const{
	stringstream toreturn;
	if(node->isLeaf()){
	    if(nodist)
		toreturn<<string(depth,'\t')<<node->getPopulation()<<endl;       
	    else
		toreturn<<string(depth,'\t')<<*node<<endl;       
	}else{
	  
	    toreturn<<printSubtree(node->getLeftChild(),  depth+1,nodist);
	    /* toreturn<<string(depth,'\t')<<*node<<endl;       	 */
	    toreturn<<printSubtree(node->getRightChild(), depth+1,nodist);
	}
	return toreturn.str();
    }

    friend ostream& operator<<(ostream& os, const Tree & tr){
	/* cerr<<"called"<<endl; */
	os<<tr.printSubtreeNewick(tr.root)<<";";
	return os;
    }

};
#endif
