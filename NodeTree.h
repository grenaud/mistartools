/*
 * NodeTree
 * Date: Mar-20-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef NodeTree_h
#define NodeTree_h

#include <stdlib.h>
#include <string>
#include <iostream>
#include <set>
#include "UnrootedNode.h"

using namespace std;

class NodeTree{
 private:
    bool leaf;
    bool hasParent;
    string population;
    double distance;//distance to the parent
    NodeTree * leftChild;
    NodeTree * rightChild;
    NodeTree * parent;

    void fakeNodeTreeConst(NodeTree * _leftChild, NodeTree * _rightChild,double _distance=0.0);
 public:
    NodeTree(string _population,double _distance=0.0);
    NodeTree(const UnrootedNode * unnode,const UnrootedNode * _parentNode,double _distance=0.0,bool forceDist=false); //we need to force the distance to take the _distance paramter for rooting the tree as we split in the middle
    NodeTree(NodeTree * _leftChild, NodeTree * _rightChild,double _distance=0.0);
    NodeTree(const NodeTree & other);
    ~NodeTree();
    NodeTree & operator= (const NodeTree & other);
    
    bool isLeaf() const;
    string getPopulation() const ;
    NodeTree * getLeftChild() const ;
    NodeTree * getRightChild() const ;
    NodeTree * getParent() const ;
    bool getHasParent() const ;
    
    void setParent(NodeTree * _parent)  ;
    set<string> getSubGroup() const;

    double getDistance() const ;
    void   setDistance(double dist);
    void printSubTrees() const ;

    friend ostream& operator<<(ostream& os, const NodeTree & nt){		
	os<<nt.population<<":"<<nt.distance;
	return os;
    }

};
#endif
