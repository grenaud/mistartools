/*
 * ComputeGMM
 * Date: Jul-04-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "ComputeGMM.h"

ComputeGMM::ComputeGMM(Rcpp::NumericMatrix M, int numberClasses){
//ComputeGMM::ComputeGMM(const arma::mat & M, int numberClasses){

    cerr<<"ComputeGMM"<<endl;
    cerr<<"Input matrix "<<M<<endl;

    // try {
    // 	int    argcf=0;
    // 	char * argvf [0];
    // 	R = new RInside(argcf,argvf);          // create an embedded R instance 

    // 	string txt = "suppressMessages(library(EMCluster))";
    // 	R->parseEvalQ(txt);              // load library, no return value

    // 	// (*R)["M"] = M;

    // 	// txt = "ret <- init.EM(M, nclass = "+stringify(numberClasses)+")";
    // 	// R->parseEval(txt);  // assign mat. M to NumericMatrix



    // } catch(std::exception& ex) {
    //     std::cerr << "ComputeGMM::Exception caught : " << ex.what() << std::endl;
    // } catch(...) {
    //     std::cerr << "ComputeGMM::Unknown exception caught" << std::endl;
    // }


}

ComputeGMM::~ComputeGMM(){

}


// arma::mat ComputeGMM::getMu(){
//     string txt = "ret$Mu;";
//     return R->parseEval(txt);  // assign mat. M to NumericMatrix
// }

// arma::mat ComputeGMM::probability(arma::mat T){
//     (*R)["T"] = T;
    
//     string txt = "e.step(as.matrix(T),emobj=ret)$Gamma;";
//     return R->parseEval(txt);  // assign mat. M to NumericMatrix
// }
