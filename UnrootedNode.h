/*
 * UnrootedNode
 * Date: Mar-20-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef UnrootedNode_h
#define UnrootedNode_h

#include <string>
#include <iostream>
#include <set>
#include <vector>

using namespace std;

class UnrootedNode{
 private:
    bool leaf;
    string population;
    /* double distance;//distance to the parent */
    /* UnrootedNode * leftChild; */
    /* UnrootedNode * rightChild; */
    /* UnrootedNode * parent; */
    vector<double> distance;
    vector<UnrootedNode *> neighbors;
    bool hasLeafWithNameIntern(const UnrootedNode * callingNode,string name) const;
    void printRecurIntern(const UnrootedNode * callingNode) const;
    void findNodeIntern(const UnrootedNode * callingNode,string name,UnrootedNode * & returnNode) const;

 public:
    UnrootedNode();
    UnrootedNode(string _population);
    //UnrootedNode(UnrootedNode * _leftChild, UnrootedNode * _rightChild,double _distance=0.0);

    UnrootedNode(const UnrootedNode & other);
    ~UnrootedNode();
    UnrootedNode & operator= (const UnrootedNode & other);
    
    bool isLeaf() const;
    int getNumberOfNeighbors() const;
    string getPopulation() const ;
    void addNeighbor(UnrootedNode * _neighbor,double _distance=0.0);
    bool hasLeafWithName(string name) const;
    void printRecur() const;
    void findNode(string name, UnrootedNode * & returnNode) const;
    const vector<UnrootedNode *> getNeighbors() const;
    const vector<double> getDistance() const;

    /* UnrootedNode * getLeftChild() const ; */
    /* UnrootedNode * getRightChild() const ; */
    /* UnrootedNode * getParent() const ; */
    /* void setParent(UnrootedNode * _parent)  ; */
    /* set<string> getSubGroup() const; */

    /* double getDistance() const ; */
    /* void   setDistance(double dist); */

    friend ostream& operator<<(ostream& os, const UnrootedNode & nt){		
	if(nt.leaf)
	    os<<nt.population;
	else
	    os<<"non-leaf with "<<nt.neighbors.size()<<" neighbors"<<endl;
	return os;
    }

};
#endif
