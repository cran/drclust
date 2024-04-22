#ifndef PREP_H
#define PREP_H

#include <RcppArmadillo.h>

arma::mat zscore(arma::mat X);
arma::mat minmax(arma::mat X);
bool checkPrep(int prep);
arma::mat preproc(arma::mat X, int prep);
bool checkStats(int stats);
#endif
