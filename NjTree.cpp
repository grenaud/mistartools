/*
 * NjTree
 * Date: Apr-17-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "NjTree.h"


//#define DEBUG

//Tree * neighborJoin(const vector< vector<double> >& d,vector<string>  names,vector<UnrootedNode *>  nodes,unsigned int r) {
Tree * neighborJoin(const mat & d,vector<string>  names,vector<UnrootedNode *>  nodes,unsigned int r) {
    // cerr<<"neighborJoin() "<<vectorToString(names)<<"\t"<<r<< endl;

    if (r <= 2) { //end of recursion
	//compute net divergence for all OTUs
	// vector<double> vec_r=vector<double>(r,0.0);

	// for (unsigned int i = 0; i < r; ++i) {
	//     double sum=0.0;
	//     for (unsigned int j = 0; j < r; ++j) {
	// 	sum+=d[i][j];
	//     }
	//     vec_r[i]=sum;
	// }
	
	double distI= d(0,1);
	// cout<<"distI "<<distI<<endl;
	// exit(1);
	//double distJ= d[0][1];

	UnrootedNode * left =nodes[0];
	UnrootedNode * right=nodes[1];
	
	// cerr<<"unrooted done 1"<<endl;
	// left->printRecur();
	// cerr<<"popll  "<<left->isLeaf()<<endl;
	// cerr<<"popl   "<<left->getPopulation()<<endl;
	// cerr<<endl;
	// cerr<<endl<<"----------------"<<endl;
	// right->printRecur();
	// cerr<<endl;
	// cerr<<"poprl  "<<right->isLeaf()<<endl;
	// cerr<<"popr   "<<right->getPopulation()<<endl;
	// cerr<<endl<<"unrooted done 2"<<endl;

	UnrootedNode * rootNode=0; 
	// cout<<"distI "<<distI<<endl;

	if(left->isLeaf() && right->isLeaf() ){//both leaves
	    cerr<<"Cannot build tree using only two clades"<<endl;
	    exit(1);
	}else{
	    if(!left->isLeaf() && right->isLeaf() ){
		// cerr<<"case1"<<endl;
		left->addNeighbor(right,distI);
		left->findNode("root",rootNode);
	    }else{
		if(left->isLeaf() && !right->isLeaf() ){
		    // cerr<<"case2"<<endl;
		    right->addNeighbor(left,distI);
		    right->findNode("root",rootNode);
		}else{
		    //doesn't matter who is added to whom
		    left->addNeighbor(right,distI);
		    left->findNode("root",rootNode);		
		}
	    }
	}


	// cerr<<"unrooted done 3"<<endl;
	// left->printRecur();
	// cerr<<endl<<"----------------"<<endl;


	
	if(rootNode==0){
	    cerr<<"Cannot find root on neither side"<<endl;
	    exit(1);	
	}

	if(rootNode->getNeighbors().size() != 1){
	    cerr<<"Root does not have 2 neighbors"<<endl;
	    exit(1);	
	}
	


	// cerr<<"unrooted done 1"<<endl;
	// left->printRecur();
	// cerr<<endl<<"----------------"<<endl;
	// right->printRecur();
	// cerr<<endl<<"unrooted done 2"<<endl;
	// cerr<<rootNode->getDistance().at(0)<<endl;
	//We root at half the distance between both left and right subtrees
	Tree * toreturn=new Tree(rootNode,
				 rootNode->getNeighbors().at(0),
				 rootNode->getDistance().at(0)/2.0);
	// cout<<"dist "<<toreturn->getRoot()->getLeftChild()->getDistance()<<endl;
	// cout<<"dist "<<toreturn->getRoot()->getRightChild()->getDistance()<<endl;

	// exit(1);
	return toreturn;
    } else {
	//compute net divergence for all OTUs
	// cout<<"pre vec_R"<<endl;

	//vector<double> vec_r=vector<double>(r,0.0);
	vec vec_r (r);
	for (unsigned int i = 0; i < r; ++i) {
	    double sum=0.0;
	    for (unsigned int j = 0; j < r; ++j) {
		sum+=d(i,j);
	    }
	    vec_r(i)=sum;
	}

	//	cout<<"post vec_R"<<endl;
	double minM=DBL_MAX;
	unsigned int minMI=0;
	unsigned int minMJ=0;

	//vector< vector<double> > m(r, vector<double>(r,0.0));
	mat m (r,r);
	for (unsigned int i = 0; i < r; ++i) {
	    
	    for (unsigned int j = 0; j < r; ++j) {
		if(i==j){
		    m(i,j) = 0.0;
		}else{
		    m(i,j) = d(i,j) - (  (vec_r(i) + vec_r(j)) / double(r-2) );
		    if(m(i,j) < minM){
			minM=m(i,j);
			minMI=i;
			minMJ=j;
		    }
		}
	    }	    
	}

	// cout<<"post minM"<<endl;
	// cout<<"r "<<r<<endl;
	// cout<<"d "<<d.size()<<endl;

	// //construct new distance matrix with size -1, remove 2 OTU, add a new one
	mat newD=d;

	newD.shed_col(max(minMI,minMJ));
	newD.shed_row(max(minMI,minMJ));

	newD.shed_col(min(minMI,minMJ));
	newD.shed_row(min(minMI,minMJ));
	

// 	vector< vector<double> > newD (r-1, vector<double>(r-1,0.0));//map_distance_matrix(d, minQA, minQB);
// 	for (unsigned int i = 0; i < r; i++) {	    
// 	    for (unsigned int j = 0; j < r; j++) {

// 		if(i==minMI || i==minMJ)
// 		    continue;
// 		if(j==minMI || j==minMJ)
// 		    continue;
// 		unsigned i_2=i;
// 		unsigned j_2=j;
// 		if(i>minMI)
// 		    i_2--;
// 		if(i>minMJ)
// 		    i_2--;
// 		if(j>minMI)
// 		    j_2--;
// 		if(j>minMJ)
// 		    j_2--;

// #ifdef DEBUG
// 		cout<<"a) "<<i<<"\t"<<j<<endl;
// 		cout<<"a) "<<i_2<<"\t"<<j_2<<endl;
// 		cout<<"a) "<<minMI<<"\t"<<minMJ<<endl;

// 		cout<<d[i][j]<<endl;
// 		cout<<"b) "<<i<<"\t"<<j<<endl;
// 		cout<<newD[i_2][j_2] <<endl;
// #endif
		
// 		newD[i_2][j_2] = d[i][j] ;
// 	    }    
// 	}



	//cout<<"-----pre new OTU-------"<<endl;
	// #ifdef DEBUG
	// 	for(unsigned int i=0;i<(r-1);i++){
	// 	    cout<<vectorToString(newD[i],"\t")<<endl;
	// 	}
	// #endif

	//new distance for new OTU
	newD.resize(newD.n_rows+1,newD.n_cols+1);
	
	unsigned int k2=0;
	for(unsigned int k=0;k<r;k++){
	    if(k==minMI || k==minMJ)
		continue;
	    double newdist=((d(minMI,k)+d(k,minMJ) - d(minMI,minMJ))/2.0);	    
	    //cout<<names[k]<<"\t"<<names[minMI]<<"\t"<<names[minMJ]<<"\t"<<newdist<<endl;
	    newD(r-2,k2) = newdist;
	    newD(k2,r-2) = newdist;
	    k2++;
	}
	newD(r-2,r-2)=0.0;


	double distI=( (d(minMI,minMJ))/2.0 + (vec_r(minMI)-vec_r(minMJ))/( 2.0*double(r-2) ) );
	double distJ=( (d(minMI,minMJ))/2.0 + (vec_r(minMJ)-vec_r(minMI))/( 2.0*double(r-2) ) );

	UnrootedNode * newnode=new UnrootedNode();
	newnode->addNeighbor( nodes[minMI] ,distI);
	newnode->addNeighbor( nodes[minMJ] ,distJ);

	// cerr<<"new node "<<endl;
	// newnode->printRecur();
	// cerr<<endl<<"new node "<<endl;

	//erasing from names
	string newname="";
	if(minMI<minMJ){
	    newname="("+names[minMI]+","+names[minMJ]+")";
	    names.erase(names.begin()+minMJ);
	    names.erase(names.begin()+minMI);
	    nodes.erase(nodes.begin()+minMJ);
	    nodes.erase(nodes.begin()+minMI);
	}else{
	    newname="("+names[minMI]+","+names[minMJ]+")";
	    names.erase(names.begin()+minMI);
	    names.erase(names.begin()+minMJ);
	    nodes.erase(nodes.begin()+minMI);
	    nodes.erase(nodes.begin()+minMJ);
	}
	names.push_back(newname);
	nodes.push_back(newnode);

	


	return neighborJoin(newD,names,nodes, r-1);
    }
    
    return 0;
}

// Tree * makeTree(vector< vector<DistanceResult> > dr,const vector<string> * names,int numberOfPopulations){
    
    
//     vector< vector<double> > distanceMatrix = vector< vector<double> > (numberOfPopulations, vector<double>(numberOfPopulations));

//     for(int i=0;i<numberOfPopulations;i++){
// 	for(int j=0;j<numberOfPopulations;j++){

// 	    if(i<j)
// 		distanceMatrix[i][j]=double(dr[i][j].all.getMutations());
// 	    else
// 		distanceMatrix[i][j]=double(dr[j][i].all.getMutations());
// 	}
//     }



//     // for(int i=0;i<numberOfPopulations;i++){
//     // 	for(int j=0;j<numberOfPopulations;j++){
//     // 	    cout<<names->at(i)<<"\t"<<names->at(j)<<"\t"<<distanceMatrix[i][j]<<endl;
//     // 	}
//     // }

//     return neighborJoinFromDist(distanceMatrix,names,numberOfPopulations);
// }

Tree * neighborJoinFromDist(const vector< vector<double> >& distanceMatrix,const vector<string> &  names,int numberOfPopulations) {
    return neighborJoinFromDist(vectorDouble2mat(distanceMatrix),names,numberOfPopulations) ;
}

Tree * neighborJoinFromDist(const mat & dr,const vector<string> &  names,int numberOfPopulations) {
    //    cout<<"neighborJoinFromDist()"<<endl;
    vector<string> namesToUse = vector<string>(names);
    vector<UnrootedNode *> nodes;
    if(dr.n_cols != names.size()){
	cerr<<"The distance matrix is not equal in size to the names"<<endl;
	exit(1);
    }

    if(dr.n_rows != names.size()){
	cerr<<"The distance matrix is not equal in size to the names"<<endl;
	exit(1);
    }


#ifdef DEBUG
	// cout<<"\t"<<vectorToString(names,"\t")<<endl;
	// for(int i=0;i<numberOfPopulations;i++){
	//     vector<double> temp;
	//     for(int j=0;j<numberOfPopulations;j++){
	// 	temp.push_back(distanceMatrix[i][j]);
	//     }
	//     cout<<names[i]<<"\t"<<vectorToString(temp,"\t")<<endl;
	// }
#endif

    // cout<<"neighborJoinFromDist2()"<<endl;

    for(int i=0;i<numberOfPopulations;i++){
	UnrootedNode * newnode=new UnrootedNode(names[i]);
	nodes.push_back(newnode);
    }
    // cout<<"neighborJoinFromDist3()"<<endl;

    
    Tree * toReturn = neighborJoin(dr,namesToUse,nodes,numberOfPopulations);
    // cout<<"neighborJoinFromDist4()"<<endl;

    return toReturn;

}
