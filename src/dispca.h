#ifndef DISPCA_H
#define DISPCA_H

#include <RcppArmadillo.h>

arma::mat dpca_updateA(int posmax, int g, arma::vec ibCg, arma::vec ibCpm, arma::mat S, arma::mat A, arma::vec JJ, int J, int n);

bool checkConstr(arma::vec constr, int J, int Q);


#endif
