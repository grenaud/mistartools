/*
 * UnrootedNode
 * Date: Mar-20-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "UnrootedNode.h"


UnrootedNode::UnrootedNode(){
    leaf=true;
}


UnrootedNode::UnrootedNode(string _population){
    leaf=true;
    population=_population;
    //distance=_distance;
}

void UnrootedNode::addNeighbor(  UnrootedNode * _neighbor,double _distance){
    leaf=false;   
    neighbors.push_back(_neighbor);
    distance.push_back(_distance);
    //adding back reference
    _neighbor->neighbors.push_back(this);
    _neighbor->distance.push_back(_distance);
}



const vector<UnrootedNode *>  UnrootedNode::getNeighbors() const{
    return neighbors;
}

const vector<double>  UnrootedNode::getDistance() const{
    return distance;
}


bool UnrootedNode::isLeaf() const{
    return leaf;
}


string UnrootedNode::getPopulation() const{
    return population;
}


bool UnrootedNode::hasLeafWithNameIntern(const UnrootedNode * callingNode,string name) const{
    //cout<<"hasLeafWithNameIntern "<<
    if(leaf){
	return (population==name);
    }else{
	bool toReturn=false;
	for(unsigned int i=0;i<neighbors.size();i++){
	    if(callingNode != neighbors[i])
		toReturn = toReturn || neighbors[i]->hasLeafWithNameIntern(this,name);
	}
	return toReturn;
    }
}

bool UnrootedNode::hasLeafWithName(string name) const{
    return hasLeafWithNameIntern(this,name);
}

int UnrootedNode::getNumberOfNeighbors() const{
    return neighbors.size();
}

void UnrootedNode::findNodeIntern(const UnrootedNode * callingNode,string name, UnrootedNode * &returnNode) const{

    if(leaf){
	// cout<<"findNodeIntern "<<population<<endl;
	// cout<<returnNode<<endl;
	// if(population==name)
	//     returnNode=this;
	// cout<<returnNode<<endl;
    }else{
	for(unsigned int i=0;i<neighbors.size();i++){
	    if(callingNode != neighbors[i]){
		if(neighbors[i]->isLeaf() && neighbors[i]->getPopulation() == name){
		    //cout<<"findNodeIntern found "<<population<<endl;
		    returnNode=neighbors[i];
		    //cout<<returnNode<<endl;
		}else{
		    neighbors[i]->findNodeIntern(this,name,returnNode);
		}
	    }
	}
    }
}

void UnrootedNode::findNode(string name, UnrootedNode * & returnNode) const{
    findNodeIntern(this,name,returnNode);
}


void UnrootedNode::printRecurIntern(const UnrootedNode * callingNode) const{
    //cout<<"hasLeafWithNameIntern "<<
    if(leaf){
	cout<<population;
    }else{
	cout<<"(";
	bool printedParant=false;
	for(unsigned int i=0;i<neighbors.size();i++){

	    if(callingNode != neighbors[i]){ //if the calling node is not this neighbor, we recurse
		neighbors[i]->printRecurIntern(this);
		cout<<":"<<distance[i];
		if(!printedParant){
		    printedParant=true;
		    cout<<",";
		}
	    }

	}
	cout<<")";
    }
}

void UnrootedNode::printRecur() const{
    return printRecurIntern(this);
}





// UnrootedNode::UnrootedNode(  UnrootedNode * _leftChild,UnrootedNode * _rightChild,double _distance){
//     leaf=false;

//     leftChild  = _leftChild;
//     rightChild = _rightChild;

//     if(leftChild->getPopulation() < rightChild->getPopulation()){
// 	population="("+leftChild->getPopulation()+ "," + rightChild->getPopulation()+")";
//     }else{
// 	population="("+rightChild->getPopulation()+ "," + leftChild->getPopulation()+")";
//     }

//     distance=_distance;
//     leftChild->setParent(this);
//     rightChild->setParent(this);

// }

// UnrootedNode::~UnrootedNode(){
//     if(leaf){
// 	delete(leftChild);
// 	delete(rightChild);    
//     }
// }





// UnrootedNode * UnrootedNode::getLeftChild() const {
//     return leftChild;
// }

// UnrootedNode * UnrootedNode::getRightChild() const {
//     return rightChild;
// }


// double UnrootedNode::getDistance() const {
//     return distance;
// }

// void UnrootedNode::setDistance(double dist){
//     distance=dist;
// }


// UnrootedNode * UnrootedNode::getParent() const {
//     return parent;
// }

// void UnrootedNode::setParent(UnrootedNode * _parent)  {
//     parent=_parent;
// }


// set<string> UnrootedNode::getSubGroup() const {
//     set<string> toreturn;

//     if(leaf){
// 	toreturn.insert(population);
//     }else{
// 	set<string> setLeft  =  leftChild->getSubGroup();
// 	set<string> setRight = rightChild->getSubGroup();
// 	toreturn.insert(setLeft.begin(),   setLeft.end());
// 	toreturn.insert(setRight.begin(), setRight.end());
//     }

//     return toreturn;
// }


