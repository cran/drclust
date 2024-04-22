#include <RcppArmadillo.h>
//using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

//' @name mrand
//' @title Adjusted Rand Index
//' @description
//' Performs the Adjusted Rand Index on a confusion matrix (row-by-column product of two partition-matrices). ARI is a measure of the similarity between two data clusterings.
//' 
//' @usage mrand(N)
//' 
//' @param N Confusion matrix.
//' 
//' @return \item{mri}{Adjusted Rand Index of a confusion matrix (scalar).}
//' 
//' @author Ionel Prunila, Maurizio Vichi
//' 
//' @references Rand W. M. (1971) "Objective criteria for the evaluation of clustering methods" <doi:10.2307/2284239>
//' 
//' @examples
//' # Iris data 
//' # Loading the numeric variables of iris data
//' iris <- as.matrix(iris[,-5]) 
//' 
//' # standardizing the data
//' iris <- scale(iris)
//' 
//' # double k-means with 3 unit-clusters and 2 components for the variables
//' p1 <- redkm(iris, K = 3, Q = 2, Rndstart = 10)
//' p2 <- doublekm(iris, K=3, Q=2, Rndstart = 10)
//' mri <- mrand(t(p1$U)%*%p2$U)
//' @export

//[[Rcpp::export]]

double mrand(arma::mat N){
  double mri;
  int n = arma::sum(arma::sum(N));
  double sumi = 0.5*(arma::sum(arma::pow(arma::sum(N.t()),2))-n);
  double sumj = 0.5*(arma::sum(arma::pow(arma::sum(N),2))-n);
  double pb = sumi*sumj/(n*(n-1)/2);
  return mri = (0.5*(arma::sum(arma::sum(arma::pow(N,2)))-n)-pb)/((sumi+sumj)/2-pb);
}

