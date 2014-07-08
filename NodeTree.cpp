/*
 * NodeTree
 * Date: Mar-20-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "NodeTree.h"

NodeTree::NodeTree(string _population,double _distance){
    hasParent=false;
    leaf=true;
    population=_population;
    distance=_distance;
}

NodeTree::NodeTree(const UnrootedNode * unnode,const UnrootedNode * _parentNode,double _distance,bool forceDist){
    hasParent=false;
    // cout<<endl<<"NodeTree constructor"<<endl;
    if(unnode->isLeaf()){
	//	cout<<endl<<"NodeTree leaf "<<unnode->getPopulation()<<endl;
	leaf=true;
	population=unnode->getPopulation();
	distance=_distance;
    }else{
	//cout<<endl<<"NodeTree recurs "<<endl;
	//unnode->printRecur();
	UnrootedNode * leftUnNode =0;
	UnrootedNode * rightUnNode=0;
	double lDist=0.0;
	double rDist=0.0;
	double distanceToParent=0;

	for(unsigned int i=0;i<unnode->getNeighbors().size();i++){
	    //setting the left node and right node
	    if(_parentNode != unnode->getNeighbors().at(i)){
		if(leftUnNode == 0){
		    leftUnNode=unnode->getNeighbors().at(i);
		    lDist     =unnode->getDistance().at(i);
		}else{
		    if(rightUnNode == 0){
			rightUnNode=unnode->getNeighbors().at(i);
			rDist      =unnode->getDistance().at(i);
		    }
		}
	    }else{ //setting the parent node
		//distanceToParent=unnode->getDistance().at(i);
		//use to root the tree in the middle
		//the length of the branch to the parent is cut in two
		//therefore no sense in taking  unnode->getDistance().at(i);
		if(forceDist){
		    distanceToParent = _distance;
		}else{ 
		    distanceToParent = unnode->getDistance().at(i);
		}
	    }	    
	}


	if(leftUnNode == 0 || rightUnNode==0){
	    cerr<<"NodeTree found non-leaf node with less than 2 children node"<<endl;
	    exit(1);
	}

	// cout<<endl<<"recurs2 left:"<<endl;
	// leftUnNode->printRecur();
	// cout<<endl<<"recurs2 right:"<<endl;
	// rightUnNode->printRecur();

	NodeTree * _leftChild  = new NodeTree(leftUnNode, unnode,lDist);
	NodeTree * _rightChild = new NodeTree(rightUnNode,unnode,rDist);
	//this=new NodeTree(_leftChild,_rightChild,distance);
	//cout<<"new dist "<<distance<<endl;
	fakeNodeTreeConst(_leftChild,_rightChild,distanceToParent);
	//cout<<"end const node1"<<endl;
	//printSubTrees();
	//cout<<"end const node2"<<endl;
    }
}

NodeTree::NodeTree( NodeTree * _leftChild,NodeTree * _rightChild,double _distance){
    hasParent=false;
    fakeNodeTreeConst( _leftChild,_rightChild,_distance);
}

void NodeTree::fakeNodeTreeConst( NodeTree * _leftChild,NodeTree * _rightChild,double _distance){
    leaf=false;

    leftChild  = _leftChild;
    rightChild = _rightChild;

    if(leftChild->getPopulation() < rightChild->getPopulation()){
	population="("+leftChild->getPopulation()+ "," + rightChild->getPopulation()+")";
    }else{
	population="("+rightChild->getPopulation()+ "," + leftChild->getPopulation()+")";
    }

    distance=_distance;
    leftChild->setParent(this);
    rightChild->setParent(this);
}

NodeTree::~NodeTree(){
    if(leaf){
	delete(leftChild);
	delete(rightChild);    
    }
}


bool NodeTree::getHasParent() const{
    return hasParent;
}

bool NodeTree::isLeaf() const{
    return leaf;
}


string NodeTree::getPopulation() const{
    return population;
}


NodeTree * NodeTree::getLeftChild() const {
    return leftChild;
}

NodeTree * NodeTree::getRightChild() const {
    return rightChild;
}


double NodeTree::getDistance() const {
    return distance;
}

void NodeTree::setDistance(double dist){
    distance=dist;
}


NodeTree * NodeTree::getParent() const {
    return parent;
}

void NodeTree::setParent(NodeTree * _parent)  {
    hasParent=true;
    parent=_parent;
}


set<string> NodeTree::getSubGroup() const {
    set<string> toreturn;

    if(leaf){
	toreturn.insert(population);
    }else{
	set<string> setLeft  =  leftChild->getSubGroup();
	set<string> setRight = rightChild->getSubGroup();
	toreturn.insert(setLeft.begin(),   setLeft.end());
	toreturn.insert(setRight.begin(), setRight.end());
    }

    return toreturn;
}

void NodeTree::printSubTrees() const {


    if(leaf){
	cout<<"LEAF:"<<population<<endl;
    }else{
	cout<<"(";
	leftChild->printSubTrees();
	cout<<",";
	rightChild->printSubTrees();
	cout<<")";
    }

}



