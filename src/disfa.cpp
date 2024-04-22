#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>
//#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


arma::vec ACP(arma::mat Xr){
  arma::vec y;
  int maxit = 300;
  int n = Xr.n_rows;
  int Q = Xr.n_cols;
  arma::vec a = arma::randu(Q,1);
  double tol = 1e-9;
  double error = arma::datum::inf;
  double last = arma::datum::inf;
  int itr = 0;
  arma::vec e;
  while(!(std::abs(last-error)<error*tol) && itr <= maxit){
    itr++;
    y = (Xr*a)*(1.0/(a.t()*a));
    a = (Xr.t()*y)*(1.0/(y.t()*y));
    last = error;
    e = y-(Xr*a);
    error = arma::conv_to<double>::from((e.t()*e)/n);
  }
  a = a*(1.0/arma::sqrt(a.t()*a));
  y = Xr*a;
  return a;
}
//
arma::vec AF1(arma::mat S){
  arma::vec U;
  int maxit = 300;
  arma::vec dS = arma::diagvec(S);
  arma::mat Psi = arma::diagmat(dS);
  
  //Given initial Psi, compute A
  arma::mat Psii = arma::diagmat(arma::pow(1.0/dS,0.5));
  arma::mat D = Psii*S*Psii;
  U = ACP(D);
  arma::mat L = U.t()*D*U;
  arma::mat Ps = arma::diagmat(pow(dS,0.5));
  arma::vec ll = arma::diagvec(L);
  ll = arma::pow((ll-1),0.5);
  arma::vec A = Ps*U*arma::diagmat(ll);
  
  //discrepancy
  arma::mat AA = A*A.t();
  arma::mat Sx = AA + Psi;
  double discrep0 = std::log(arma::det(Sx)) + arma::trace(/*arma::pinv(Sx)*S*/solve(Sx,S, arma::solve_opts::allow_ugly));
  
  double tol = 1e-10;
  double discrep = arma::datum::inf;
  int itr = 0;
  arma::vec dAA, dPsi;
  while(!(std::abs(discrep - discrep0)<tol) && itr<=maxit){
    itr++;
    if(itr!= 1)
      discrep0 = discrep;
    // given A, update Psi
    AA = A*A.t();
    dAA = arma::diagvec(AA);
    Psi = arma::diagmat(dS-dAA);
    
    // given Psi, update A
    dPsi = arma::diagvec(Psi);
    Psii = arma::diagmat(arma::pow((1.0/dPsi),0.5));
    D = Psii*S*Psii;
    U = ACP(D);
    L = U.t()*D*U;
    Ps = arma::diagmat(arma::pow(dPsi,0.5));
    ll = arma::diagvec(L);
    ll = ll-1;
    ll = arma::pow(ll,0.5);
    A = Ps*U*arma::diagmat(ll);
    AA = A*A.t();
    Sx = AA + Psi;
    discrep = std::log(arma::det(Sx)) + arma::trace(/*arma::pinv(Sx)*S*/solve(Sx,S, arma::solve_opts::allow_ugly));
  }
  return A;
}


//' @name CronbachAlpha
//' @title Cronbach Alpha
//' @description
//' Computes the Cronbach Alpha index on a units x variables data matrix. It measures the internal reliability, i.e., the propensity of J variables of a data matrix (n units x J variables) to be concordantly correlated with a single factor (composite indicator).
//' 
//' @usage CronbachAlpha(X)
//' 
//' @param X Units x variables numeric data matrix.
//' 
//' @return \item{as}{Cronbach's Alpha}
//' 
//' @author Ionel Prunila, Maurizio Vichi
//' 
//' @references Cronbach L. J. (1951) "Coefficient alpha and the internal structure of tests" <doi:10.1007/BF02310555>
//' 
//' @examples
//' # Iris data 
//' # Loading the numeric variables of iris data
//' iris <- as.matrix(iris[,-5]) 
//' 
//' # standardizing the data
//' iris <- scale(iris)
//' 
//' # compute Cronbach's Alpha
//' as <- CronbachAlpha(iris)
//' @export

//[[Rcpp::export]]
double CronbachAlpha(arma::mat X){
  int nItem = X.n_cols;
  arma::mat B = arma::trimatu(arma::ones(nItem,nItem),1);
  arma::mat r = arma::cor(X);
  double m_r = arma::mean(r.elem(arma::find(B==1)));
  double as;
  as = nItem*m_r/(1+((nItem-1)*m_r));
  return as;
}

double r2pv(double r, int n){
  double p;
  if(n<3)
    stop("Error: n < 3");
  if(r==1)
    p=0;
  double t = std::sqrt(n-2)*r/(std::sqrt(1-r*r));
  t = std::abs(t);
  p = 2*(1 - R::pt(t, n-2, true, false));
  return p;
}

