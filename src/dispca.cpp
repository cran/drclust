#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>
//#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat dpca_updateA(int posmax, int g, arma::vec ibCg, arma::vec ibCpm, arma::mat S, arma::mat A, arma::vec JJ, int J, int n){
  
  arma::vec l, ll, zJ = arma::zeros(J,1);
  arma::mat Av, AAv;
  int a, aa;
  
  arma::uvec indices = arma::find(ibCg==1);
  arma::vec JCg = JJ.elem(indices);
  
  indices = arma::find(ibCpm==1);
  arma::vec JCpm = JJ.elem(indices);
  
  arma::mat Sg = S.submat(arma::conv_to<arma::uvec>::from(JCg), arma::conv_to<arma::uvec>::from(JCg));
  arma::mat S0g = S.submat(arma::conv_to<arma::uvec>::from(JCpm), arma::conv_to<arma::uvec>::from(JCpm));
  
  if(sum(ibCg)>1){
    if(n > 20000)
      arma::eig_sym(l, Av, Sg, "dc");
    else
      arma::eig_sym(l, Av, Sg);
    a = arma::index_max(abs(l));
    A.col(g) = zJ;
    if(arma::sum(Av.col(a))<0){
      Av.col(a) = -Av.col(a);
    }
    A.submat(arma::conv_to<arma::uvec>::from(JCg), arma::uvec{static_cast<unsigned int>(g)}) = Av.col(a);
  }
  else{
    A.col(g) = zJ;
    A.submat(arma::conv_to<arma::uvec>::from(JCg), arma::uvec{static_cast<unsigned int>(g)}) = arma::ones(JCg.n_elem);
  }
  if(arma::sum(ibCpm)>1){
    if(n > 20000)
      arma::eig_sym(ll, AAv, S0g, "dc");
    else
      arma::eig_sym(ll, AAv, S0g);
    aa = arma::index_max(abs(ll));
    if(arma::sum(AAv.col(aa))<0){
      AAv.col(aa) = -AAv.col(aa);
    }
    A.col(posmax) = zJ;
    A.submat(arma::conv_to<arma::uvec>::from(JCpm), arma::uvec{static_cast<unsigned int>(posmax)}) = AAv.col(aa);
  }
  else{
    A.col(posmax) = zJ;
    A.submat(arma::conv_to<arma::uvec>::from(JCpm), arma::uvec{static_cast<unsigned int>(posmax)}) = arma::ones(JCpm.n_elem);
  }
  return A;
}

bool checkConstr(arma::vec constr, int J, int Q){
  bool args_ok = 1;
  int c = constr.size();
  if(c/*onstr.size()*/<J && constr.size()!=1){
    args_ok = 0;
    warning("Error: the length of the constr vector must be = J = nr. of columns");
  }
  else{
    for(int j=0;j<c/*onstr.size()*/;j++){
      if(constr(j)<0 || constr(j)>Q){
        args_ok = 0;
        warning("Error: values of 'constr' vector should be >=1 and <= Q");
        break;
      }
    }
  }
  return args_ok;
}
