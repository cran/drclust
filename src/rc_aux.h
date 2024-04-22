#ifndef RC_AUX_H
#define RC_AUX_H

#include <RcppArmadillo.h>

arma::mat randPU(int n, int c);
bool checkArgs(int Q, int Rndstart, int stats, int maxiter, double eps, int J);
bool checkK(int K, int n);
arma::mat km(arma::mat X, int K, int Rndstart);
arma::mat split_maxwd(arma::rowvec su, arma::rowvec wd, arma::mat U, arma::mat Xs, int K);
arma::mat assign_ssed(arma::mat S2x, arma::mat M, arma::mat Xs, int K, int n);

#endif
