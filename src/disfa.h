#ifndef DISFA_H
#define DISFA_H

#include <RcppArmadillo.h>

arma::vec ACP(arma::mat Xr);
arma::vec AF1(arma::mat S);
double CronbachAlpha(arma::mat X);
double r2pv(double r, int n);
#endif
