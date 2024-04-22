#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

//' @name cluster
//' @title classification variable
//' @description
//' Recodes the binary and row-stochastic membership matrix U into the classification variable (similar to the "cluster" output returned by kmeans()).
//' 
//' @usage cluster(U)
//' 
//' @param U Binary and row-stochastic matrix.
//' 
//' @return \item{cl}{vector of length n indicating, for each element, the index of the cluster to which it has been assigned.}
//' 
//' @author Ionel Prunila, Maurizio Vichi
//' 
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
//' p1 <- redkm(iris, K = 3, Q = 2)
//' cl <- cluster(p1$U)
//' 
//' @export

//[[Rcpp::export]]

arma::rowvec cluster(arma::mat U){
  unsigned int n = U.n_rows; 
  arma::rowvec clust(n);
  for(unsigned int i = 0; i<n; i++){
    clust(i) = arma::index_max(U.row(i))+1;
  }
  return clust;
}

