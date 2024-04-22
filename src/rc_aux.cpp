#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat randPU(int n, int c) {
  arma::mat U(n, c, arma::fill::zeros);
  U.rows(0,c-1) = arma::eye(c,c);
  U(arma::span(c,n-1), 0) = arma::ones(n-c);
  for (int i = c; i < n; i++) {
    U.row(i) = arma::shuffle(U.row(i));
  }
  U = arma::shuffle(U, 0);
  return U;
}

bool checkArgs(int Q, int Rndstart, int verbose, int maxiter, double eps, int J){
  bool args_ok = 1; 
  if(Q<1 || Q>J-1){ //Col Clusters
    args_ok = 0;
    warning("Q must be specified as an integer > 0 and < J = nr. of variables");
  }
  if(Rndstart<1 || Rndstart>1000){ //Random Starts
    args_ok = 0;
    warning("Rndstart must be an integer > 0 and < 1000");
  }
  if(eps>=0.1 || eps<=0){ //Error Tollerance
    args_ok = 0;
    warning("eps must a value > 0 and < 0.1");
  }
  if(verbose!=0 && verbose != 1){ //Stat prints
    args_ok = 0;
    warning("verbose must be a value either = 0 or = 1");
    
  }
  if(maxiter < 0 || maxiter > 1000){ // Nr. of iterations
    args_ok = 0;
    warning("maxiter must be an integer > 0 and < 1000");
  }
  
  return args_ok;
}

bool checkK(int K, int n){
  bool args_ok = 1;
  if(K<1 || K>=n){ // Row Clusters
    args_ok = 0;
    warning("K must be an integer => 1 and < n = nr. of observations");
  }
  return args_ok;
}

arma::mat km(arma::mat X, int K, int Rndstart) {
  int n = X.n_rows;
  int J = X.n_cols;

  arma::mat S2x = arma::repmat(arma::sum(arma::pow(X,2),1),1,K);
  arma::mat U = randPU(n,K);
  arma::rowvec su = arma::sum(U,0);
  arma::mat M = arma::repmat((1.0/su).t(),1,J)%(U.t()*X);
  arma::vec s2m = arma::sum(arma::pow(M,2),1);

  double ssbold = arma::sum(s2m.t()%su);
  int it = 0;
  double dif = arma::datum::inf;
  double ssb;
  arma::mat Dist;
  int sumin, sumax, o;
  while(dif > 1e-8){
    it++;
    Dist = S2x + arma::repmat(s2m.t(),n,1) - 2*X*M.t();
    U = arma::conv_to<arma::mat>::from(arma::repmat(arma::min(Dist,1),1,K) == Dist);
    while(arma::sum(arma::sum(U,1)>1)>0){//check for not-unique assigments
    	arma::vec surows = arma::sum(U,1); 
    	sumax = arma::index_max(surows>1);
    	int posmax = arma::index_max(U.row(sumax));
    	U.row(sumax).fill(0);
    	U.submat(arma::uvec{static_cast<unsigned int>(sumax)}, arma::uvec{static_cast<unsigned int>(posmax)}).fill(1);
  }
    su = arma::sum(U,0);
    while(arma::sum(su==0)>0){
      sumin = arma::index_min(su);
      sumax =  arma::index_max(su);
      o = arma::index_max(U.col(sumax));
      U(o, sumin) = 1;
      U(o, sumax) = 0;
      su = arma::sum(U,0);
    }

    //centroids
    M = arma::repmat((1.0/su).t(),1,J)%(U.t()*X);
    s2m = arma::sum(arma::pow(M,2),1);

    //check convergence
    ssb = arma::sum(s2m.t()%su);
    dif = ssb - ssbold;
    ssbold = ssb;
  }
  return U;
}


arma::mat split_maxwd(arma::rowvec su, arma::rowvec wd, arma::mat U, arma::mat Xs, int K){
  unsigned int sumin = arma::index_min(su);
  unsigned int wdmax =  arma::index_max(wd);
  arma::uvec splitk = {sumin, wdmax};
  arma::uvec isplit = arma::find(U.col(wdmax)==1);
  U.rows(isplit) = arma::zeros(isplit.n_elem,K);

  if(isplit.n_elem >2){
    U.submat(isplit, splitk) = km(Xs.rows(isplit),2,10);
  }
  else if(isplit.n_elem == 2){
    U.submat(isplit, splitk) = arma::eye(2,2);
  }
  return U;
}


arma::mat assign_ssed(arma::mat S2x, arma::mat M, arma::mat Xs, int K, int n){
  arma::vec s2m = arma::sum(arma::pow(M,2),1);
  arma::mat Dist = S2x + arma::repmat(s2m.t(),n,1) - 2*Xs*M.t();
  arma::mat U = arma::conv_to<arma::mat>::from(arma::repmat(arma::min(Dist,1),1,K) == Dist);
  arma::vec surows;
  int sumax, posmax;
  while(arma::sum(arma::sum(U,1)>1)>0){//check for not-unique assigments
    surows = arma::sum(U,1); //
    sumax = arma::index_max(surows>1);
    posmax = arma::index_max(U.row(sumax));
    U.row(sumax).fill(0);
    U.submat(arma::uvec{static_cast<unsigned int>(sumax)}, arma::uvec{static_cast<unsigned int>(posmax)}).fill(1);
  }
  return U;
}
