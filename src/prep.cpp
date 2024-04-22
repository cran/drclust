#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

arma::mat zscore(arma::mat X){
  int n = X.n_rows;
  arma::mat S = cov(X,1);
  arma::mat Xs = (X - arma::repmat(arma::mean(X,0),n,1)) * arma::pow(arma::diagmat(arma::diagvec(1/S)),0.5);
  return Xs;
}
arma::mat minmax(arma::mat X){
  int n = X.n_rows;
  arma::rowvec minX = min(X);
  arma::rowvec maxX = max(X);
  arma::mat Xs = (X - repmat(minX,n,1))/(repmat(maxX-minX, n,1));
  return Xs;
}
arma::mat preproc(arma::mat X, int prep){
  arma::mat Xs;
  if(prep == 0)
    Xs = X;
  else if(prep == 1)
      Xs = zscore(X);
  else if(prep == 2)
    Xs = minmax(X);
  return Xs;
}

bool checkPrep(int prep){
  bool args_ok = 1;
  if(prep != 1 && prep!= 2 && prep != 0){
    Rcpp::Rcout << "Prep must be either 0 (no pre-processing), 1 (z-transform) or 2 (minmax)." << std::endl;
    args_ok = 0;
  }
  return args_ok;
}

bool checkStats(int stats){
  bool args_ok = 1;
  if(stats!= 0 && stats != 1){
    Rcpp::Rcout << "stats must be either 0 or 1." << std::endl;
    args_ok = 0;
  }
  return args_ok;
}
