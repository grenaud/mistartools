/*
 * mistarArmadillo
 * Date: Apr-01-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef mistarArmadillo_h
#define mistarArmadillo_h

#include <vector>
/* #include <mlpack/core.hpp> */
#include <armadillo>
/* #include <RInside.h>                  */

#include "MistarParser.h"
#include "GenomicRange.h"

using namespace arma;
using namespace std;

mat vectorDouble2mat(const vector< vector<double>  > & inputVec);
/* Rcpp::NumericMatrix mat2Rmat(const arma::mat & M); */

#endif
