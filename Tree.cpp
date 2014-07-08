/*
 * Tree
 * Date: Mar-21-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "Tree.h"

Tree::Tree(NodeTree * _root){
    root = _root;
    allSubGroups = new vector< set<string> >  ( returnAllSubgroup() );
}

Tree::Tree(const UnrootedNode * _root,const UnrootedNode * _rest,double distance){

    // cout<<"tree constructor11 "<<endl;
    // cout<<distance/2.0<<endl;
    NodeTree * lefttree  = new NodeTree(_root,_rest,distance,true); 
    // cout<<lefttree->getDistance()<<endl;
    // cout<<"tree constructor12 "<<endl;
    // cout<<distance/2.0<<endl;
    NodeTree * righttree = new NodeTree(_rest,_root,distance,true);
    // cout<<righttree->getDistance()<<endl;
    // cout<<"tree constructor "<<endl;
    // _rest->printRecur();
    // cout<<"tree constructor2 "<<endl;
    root = new NodeTree(lefttree,righttree,0);
    // cout<<"tree constructor3 "<<endl;
    // cout<<lefttree<<endl;
    // cout<<righttree<<endl;

    // cout<<"tree constructor4"<<endl;
    // cout<<*this<<endl;
    // cout<<"tree constructor5"<<endl;


}


Tree::~Tree(){
    delete(root);
    delete(allSubGroups);
}


//A subtree has an order, a subgroup does not

bool Tree::hasSubTree(string treeStr){
    Tree  * tr =new Tree( parseNewick(treeStr) ); 

    bool toreturn = hasSubTreeWithName(root,tr->nameRoot());
    delete(tr);
      
    return toreturn;
}

bool Tree::hasSubTreeWithName(NodeTree * node,string subTree){

    if(node->isLeaf()){
	return false;
    }else{
	if(node->getPopulation() == subTree)
	    return true;
	else
	    return hasSubTreeWithName(node->getLeftChild(),subTree) || hasSubTreeWithName(node->getRightChild(),subTree);
    }
    
}


bool Tree::hasSubGroup(set<string> tocheck){
    for(unsigned int i=0;i<allSubGroups->size();i++){
	if(allSubGroups->at(i) == tocheck)
	    return true;
    }
    return false;
}

vector< set<string> > Tree::returnAllSubgroup(){
    vector< set<string> > toReturn;

    returnAllSubgroupNode(root,&toReturn);

    return toReturn;
}

void Tree::returnAllSubgroupNode(const NodeTree * node,vector< set<string> > * accumulateSets){
    if( node->isLeaf() ){
	accumulateSets->push_back( node->getSubGroup()  );
    }else{
	
	returnAllSubgroupNode(node->getLeftChild(),  accumulateSets);
	accumulateSets->push_back( node->getSubGroup()  );
	returnAllSubgroupNode(node->getRightChild(), accumulateSets);
    }
}


bool Tree::hasLeafWithName(string nameToLookFor){
    return hasLeafWithNameIntern(root,nameToLookFor);
}

bool Tree::hasLeafWithNameIntern(NodeTree * node,string nameToLookFor){
    if( node->isLeaf() ){
	return (node->getPopulation() == nameToLookFor);	    
    }else{
	return (hasLeafWithNameIntern(node->getLeftChild(),  nameToLookFor)|| hasLeafWithNameIntern(node->getRightChild(), nameToLookFor));
    }    
}


//goes up in tree
void Tree::getDistanceUpperTree(const NodeTree * node,
				const NodeTree * callingNode,
				vector<double> * distForNodes,
				vector<string> * distNames,
				double currentDist) const{

    cout<<"left "<<"\t"<<endl;
    node->printSubTrees();
    cout<<endl<<"left2 "<<endl;
    callingNode->printSubTrees();
    if(node->getLeftChild() != callingNode){
	cout<<"leftc "<<endl;
	node->getLeftChild()->printSubTrees();
	cout<<endl;
	getDistanceToLeave(node->getLeftChild(),
			   distForNodes,
			   distNames,
			   currentDist+node->getLeftChild()->getDistance() );
    }
    cout<<"right "<<endl;
    if(node->getRightChild() != callingNode){
	cout<<"rightc "<<endl;
	getDistanceToLeave(node->getRightChild(),
			   distForNodes,
			   distNames,
			   currentDist+node->getRightChild()->getDistance() );
    }
    cout<<"up "<<endl;
    if(node != root &&
       node->getHasParent()){
	cout<<"upc "<<endl;
	getDistanceUpperTree(node->getParent(),
			     node,			     
			     distForNodes,
			     distNames,
			     currentDist + node->getDistance() );	
    }
    cout<<"done "<<endl;
}


//goes down in tree
void Tree::getDistanceToLeave(const NodeTree * node,
			      vector<double> * distForNodes,
			      vector<string> * distNames,
			      double currentDist) const{
    
   if(node->isLeaf()){

       distNames->push_back(node->getPopulation());
       distForNodes->push_back(currentDist);
	
    }else{ 

       
       getDistanceToLeave(node->getLeftChild(),
			  distForNodes,
			  distNames,
			  currentDist+node->getLeftChild()->getDistance() );

       getDistanceToLeave(node->getRightChild(),
			  distForNodes,
			  distNames,
			  currentDist+node->getRightChild()->getDistance());
       /* toreturn<<string(depth,'\t')<<*node<<endl;       	 */
       //findInternalBranch(node->getRightChild());
    }
}



void Tree::findInternalBranch(const NodeTree * node) const{

    stringstream toreturn;
    if(node->isLeaf()){
	// if(nodist)
	//     toreturn<<string(depth,'\t')<<node->getPopulation()<<endl;
	// else
	//     toreturn<<string(depth,'\t')<<*node<<endl;
	
    }else{ 
	//if(node != root){ //we have a branch 
	vector<double>  distForNodes;
	vector<string>  distNames;

	getDistanceToLeave(node,&distForNodes,&distNames,node->getDistance()/2.0);	    
	cout<<"done getDistanceToLeave "<<endl;
	node->printSubTrees();
	cout<<endl;
	// distForNodes.clear();
	// distNames.clear();
	// getDistanceToLeave(node->getRightChild(),&distForNodes,&distNames,0.0);
	cout<<"dist "<<vectorToString(distForNodes)<<endl;
	cout<<"name "<<vectorToString(distNames)<<endl;
	    
	if(node->getHasParent()){
	    cout<<"hasParent "<<endl;
	    getDistanceUpperTree(node->getParent(),
				 node,			     
				 &distForNodes,
				 &distNames,
				 (0.0 + node->getDistance()/2.0) );
	    cout<<"hasParent done "<<endl;
	}
	cout<<"dist "<<vectorToString(distForNodes)<<endl;
	cout<<"name "<<vectorToString(distNames)<<endl;


	//}
	
	findInternalBranch(node->getLeftChild());
	/* toreturn<<string(depth,'\t')<<*node<<endl; */
	findInternalBranch(node->getRightChild());

    }

    exit(1);
    //return toreturn.str();
}


void Tree::addNewNodeUsingDist(string newName,vector<double> distForNode,vector<string> distNames){
    //for each branch, compute dist to others and insert in node with less square 
    cerr<<"addNewNodeUsingDist"<<endl;
    findInternalBranch(root->getRightChild()->getRightChild());
}

const NodeTree * Tree::getRoot() const{
    return root;
}
