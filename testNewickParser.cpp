/*
 * testNewickParser
 * Date: Mar-20-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>
#include "NewickParser.h"
#include "NodeTree.h"
#include "Tree.h"

#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {
    

    //    cout<<argv[1]<<endl;
    NewickParser np;
    Tree * tree=np.parseFile(argv[1]);
    // cout<<*tree<<endl;

    //cout<<tree->hasSubTree("(Denisova:0,DenisovaMolar1:0):0")<<"\t"<<tree->hasSubTree("((Denisova:0,DenisovaMolar1:0):0,AltNean:0):0")<<endl;
    vector< set<string> > t = tree->returnAllSubgroup();

    for(unsigned int i=0;i<t.size();i++){
	//if(t[i].size() > 1){
	    //set<string>::const_iterator it = t[i].begin();
	    //cout<<<<endl;
	    cout<<tree->hasSubGroup( t[i] )<<"\t"<<iteratorToString(t[i])<<endl;

	    //cout << *it << ",";
	    //}
	// cout<<endl;
	// }
    }
    
    // set<string> sanityCheck;
    // sanityCheck.insert("Loschbour");
    // sanityCheck.insert("Dai_B");
    // cout<<tree->hasSubGroup( sanityCheck )<<"\t"<<iteratorToString(sanityCheck)<<endl;	


    delete(tree);
    return 0;
}

