/*
 * ComputeGMM
 * Date: Jul-04-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ComputeGMM_h
#define ComputeGMM_h

/* #include <armadillo> */
/* #include <RcppArmadillo.h> */
#include <RInside.h>                    // for the embedded R via RInside


#include "utils.h" 

using namespace std;

class ComputeGMM{
private:
    Rcpp::NumericMatrix mu;
    RInside * R;
    
public:
    //ComputeGMM();
    ComputeGMM(Rcpp::NumericMatrix M, int numberClasses);
    //ComputeGMM(const arma::mat &  M, int numberClasses);

    ComputeGMM(const ComputeGMM & other);
    ~ComputeGMM();
    ComputeGMM & operator= (const ComputeGMM & other);

    /* Rcpp::NumericMatrix getMu(); */
    /* Rcpp::NumericMatrix probability(Rcpp::NumericMatrix T); */
    /* arma::mat getMu(); */
    /* arma::mat probability(arma::mat T); */



};
#endif
