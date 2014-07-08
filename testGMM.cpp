/*
 * testGMM
 * Date: Jul-05-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <armadillo>
#include <RcppArmadillo.h>
#include <RInside.h>                    // for the embedded R via RInside


#include "mistarArmadillo.h"
#include "ComputeGMM.h"
#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {
    RInside R(argc, argv);		// create an embedded R instance

    arma::mat testmat;

    testmat.load(string(argv[1]));
    cout<<testmat<<endl;

    try {
	string txt = "suppressMessages(library(EMCluster))";
	R.parseEvalQ(txt);              // load library, no return value

	R["M"] = testmat;

	txt = "ret <- init.EM(M, nclass = "+stringify(2)+")";
	R.parseEval(txt);  // assign mat. M to NumericMatrix

	txt = "ret$Mu;";
	arma::mat mu= R.parseEval(txt);
	cout<<mu<<endl ;

 	txt = "ret$pi;";
	arma::vec pi= R.parseEval(txt);
	cout<<pi<<endl ;

 	txt = "e.step(as.matrix(M),emobj=ret)$Gamma;";
	//cout<<R.parseEval(txt)<<endl;
	arma::mat P=  R.parseEval(txt); 
	cout<<P<<endl ;

    } catch(std::exception& ex) {
        std::cerr << "ComputeGMM::Exception caught : " << ex.what() << std::endl;
    } catch(...) {
        std::cerr << "ComputeGMM::Unknown exception caught" << std::endl;
    }

    // Rcpp::NumericMatrix rmat = mat2Rmat(testmat);

    // ComputeGMM * cgmm = new ComputeGMM (rmat, 2);
    // cerr<<"MU "<<endl<<cgmm->getMu()<<endl;

    // string line;
    // ifstream myFile;
    // string filename = string(argv[1]);
    // myFile.open(filename.c_str(), ios::in);

    // if (myFile.is_open()){
    //   while ( getline (myFile,line)){

    //   }
    //   myFile.close();
    // }else{
    //     cerr << "Unable to open file "<<filename<<endl;
    //     return 1;
    //  }


    return 0;
}

