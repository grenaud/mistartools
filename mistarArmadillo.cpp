/*
 * mistarArmadillo
 * Date: Apr-01-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "mistarArmadillo.h"






mat vectorDouble2mat(const vector< vector<double>  > & inputVec){
    if(inputVec.empty()){
	mat a (0,0);
	return a;
    }

    mat matToReturn (inputVec.at(0).size(),inputVec.size());
    
    for(unsigned int i=0;i<inputVec.size();i++){ 	
	vec toadd = vec( inputVec.at(i).size() );
	for(unsigned int j=0;j<inputVec.at(i).size();j++)
	    toadd[j] = inputVec.at(i).at(j);
	matToReturn.col(i)=toadd;
    }

    return matToReturn;
}

// Rcpp::NumericMatrix mat2Rmat(const arma::mat & M){
    
//     Rcpp::NumericMatrix toreturn (M.n_rows,M.n_cols);
    
//     for(unsigned int i=0;i<M.n_rows;i++){ 	

// 	for(unsigned int j=0;j<M.n_cols;j++){
// 	    toreturn(i,j) = M(i,j);
// 	}
//     }
    
//     return toreturn;
    
// }
