/*
 * NjTree
 * Date: Apr-17-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef NjTree_h
#define NjTree_h

#include <cfloat>
/* #include <mlpack/core.hpp> */
#include <armadillo>
#include "list"
#include "vector"
#include "stdlib.h"

#include "DistanceResult.h"
#include "Tree.h"
#include "NodeTree.h"
#include "UnrootedNode.h"
#include "mistarArmadillo.h"

using namespace arma;
using namespace std;

//Tree * makeTree(vector< vector<DistanceResult> > dr,const vector<string> * names,int numberOfPopulations);
Tree * makeTree(const mat & dr,const vector<string> * names,int numberOfPopulations);
/* Tree * neighborJoin(const vector< vector<double> >& d,vector<string>  names,vector<UnrootedNode *>  nodes, unsigned int r) ; */
Tree * neighborJoinFromDist(const vector< vector<double> >& distanceMatrix,const vector<string> &  names,int numberOfPopulations) ;
Tree * neighborJoinFromDist(const mat & dr,const vector<string> &  names,int numberOfPopulations) ;

/* class NjTree{ */
/* private: */

/* public: */
/* NjTree(); */
/* NjTree(const NjTree & other); */
/* ~NjTree(); */
/* NjTree & operator= (const NjTree & other); */

/* }; */
#endif
