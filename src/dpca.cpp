#define ARMA_WARN_LEVEL 0
#include "rc_aux.h"
#include "dispca.h"
#include "disfa.h"
#include "prep.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

//' @name dispca
//' @title Disjoint Principal Components Analysis
//' @description
//' Performs disjoint PCA, that is, a simplified version of PCA. Computes each one of the Q principal components from a different subset of the J variables (resulting thus, in a simplified, easier to interpret loading matrix A). 
//' 
//' 
//' @usage dispca(X, Q, Rndstart, verbose, maxiter, tol, prep, print, constr)
//' 
//' @param X Units x variables numeric data matrix.
//' @param Q Number of factors.
//' @param Rndstart Number of runs to be performed (Defaults is 20).
//' @param verbose Outputs basic summary statistics for each run (1 = enabled; 0 = disabled, default option).
//' @param maxiter Maximum number of iterations allowed (if convergence is not yet reached. Default is 100).
//' @param tol Tolerance threshold (maximum difference between the values of the objective function of two consecutive iterations such that convergence is assumed). Default is 1e-6.
//' @param prep Pre-processing of the data. 1 performs the z-score transform (default choice); 2 performs the min-max transform; 0 leaves the data un-pre-processed.
//' @param print Prints summary statistics of the results (1 = enabled; 0 = disabled, default option).
//' @param constr is a vector of length J = nr. of variables, pre-specifying to which cluster some of the variables must be assigned. Each component of the vector can assume integer values from 1 o Q (See example for more details), or 0 if no constraint on the variable is imposed (i.e., it will be assigned based on the plain algorithm).
//' 
//' @return returns a list of estimates and some descriptive quantities of the final results.
//' @return \item{V}{Variables x factors membership matrix (binary and row-stochastic). Each row is a dummy variable indicating to which cluster it has been assigned.}
//' @return \item{A}{Variables x components loading matrix.}
//' @return \item{betweenss}{Amount of deviance captured by the model (scalar).}
//' @return \item{totss}{total amount of deviance (scalar).} 
//' @return \item{size}{Number of variables assigned to each column-cluster (vector).}
//' @return \item{loop}{The index of the (best) run from which the results have been chosen.}
//' @return \item{it}{the number of iterations performed during the (best) run.}
//' 
//' @author Ionel Prunila, Maurizio Vichi
//' 
//' @references 
//' Vichi M., Saporta G. (2009) "Clustering and disjoint principal component analysis" <doi:10.1016/j.csda.2008.05.028>
//' 
//' @examples
//' # Iris data 
//' # Loading the numeric variables of iris data
//' iris <- as.matrix(iris[,-5]) 
//' 
//' # No constraint on variables
//' out <- dispca(iris, Q = 2)
//' 
//' # Constraint: the first two variables must contribute to the same factor.
//' outc <- dispca(iris, Q = 2, constr = c(1,1,0,0))
//' @export

//[[Rcpp::export]]

List dispca(Rcpp::NumericMatrix X, int Q, int Rndstart = 20, int verbose = 0, int maxiter = 100, double tol = 1e-6, int prep = 1, int print = 0, arma::vec constr = 00){
  Rcpp::List Nomi = X.attr("dimnames");
  arma::mat Xs = as<arma::mat>(X);
  /// Declaration of variables
  int J = Xs.n_cols;
  int n = Xs.n_rows;
  arma::mat Xs2, S,V, Xmean, XX;
  arma::mat L, A, A0, Y, Ymean, RotA, Sg, Av, VC = arma::eye(Q,Q), S0g, AAv, U;
  int it, loop, j, g, a, posmax, flg;
  double f0, fdif, f=0, fmax, st, varcum=0;
  arma::vec un, l, dL, varY, JJ, zJ, ibCg, JCg, ibCpm, JCpm;
  arma::uvec idL, iicc, indices, ic;
  
  // best estimates
  arma::mat Vbest, Abest, Ymbest, Ybest;
  double fbest=0, fdifbest=0;
  int loopbest, itbest;
  
  bool prep_ok = checkPrep(prep);
  Xs = preproc(Xs, prep);
  Xs2 = arma::pow(Xs,2);
  st = arma::sum(arma::sum(Xs2));
  
  S = arma::cov(Xs,1);
  bool args_ok = checkArgs(Q, Rndstart, verbose, maxiter, tol, J);
  bool constr_ok = checkConstr(constr, J, Q);
  
  
  bool stats_ok = checkStats(print);
  
  if(constr.size()==1)
    constr = arma::zeros(J,1);
  JJ = arma::regspace(0,1,J-1);
  zJ = arma::zeros(J,1);
  if(args_ok == 1 && constr_ok == 1 && prep_ok ==1 && stats_ok ==1){
    for(loop = 0; loop < Rndstart; loop++){
      it = 0;
      flg = 1;
      while(flg>0){
        V = randPU(J,Q);
        for(j=0;j<J;j++){
          if(constr(j)>0){
            V.row(j) = VC.row(constr(j)-1);
          }
        }
        flg = arma::sum(arma::find(arma::sum(V)==0));
      }
      A = arma::zeros(J,Q);
      for(g = 0; g < Q; g++){
        ibCg = V.col(g);
        JCg = JJ.elem(arma::find(ibCg==1));
        Sg = S.submat(arma::conv_to<arma::uvec>::from(JCg), arma::conv_to<arma::uvec>::from(JCg));
        if(sum(ibCg) > 1){
          if(n > 20000)
            arma::eig_sym(l, Av, Sg, "dc");
          else
            arma::eig_sym(l, Av, Sg);
          a = arma::index_max(abs(l));
          A.submat(arma::conv_to<arma::uvec>::from(JCg),arma::uvec{static_cast<unsigned int>(g)}) = Av.col(a);
        }
        else{
          A.submat(arma::conv_to<arma::uvec>::from(JCg), arma::uvec{static_cast<unsigned int>(g)}) = arma::ones(JCg.n_elem);
        }
      }
      A0=A;
      Y=Xs*A;
      f0 = arma::trace(Y.t()*Y);
      fmax = 0;
      //iteration phase
      fdif = 2*tol;
      while(it < maxiter){
        it++;
        //Given U and Ymean, update V and A
        for(j = 0; j < J; j++){
          posmax = JJ(arma::index_max(V.row(j)==1));
          if(constr(j)==0){
            for(g = 0; g < Q; g++){
              V.row(j) = VC.row(g);
              if(arma::sum(V.col(posmax))>0){
                ibCg = V.col(g); // new class of V
                ibCpm = V.col(posmax); // old class of V
                A = dpca_updateA(posmax, g, ibCg, ibCpm, S, A, JJ, J, n);
                Y = Xs*A;
                f = arma::trace(Y.t()*Y);
                if(f > fmax){
                  fmax = f;
                  posmax = g;
                  A0=A;
                }
                else{
                  A=A0;
                }
              }
            }
          }
          V.row(j) = VC.row(posmax);
        }
        Y = Xs*A;
        f = arma::trace(Y.t()*Y);
        fdif = f - f0;
        if(fdif > tol){
          f0 = f;
          fmax = f0;
          A0=A;
        }
        else{
          break;
        }
      }
      if(verbose==1)
        Rcpp::Rcout << "DPCA: Loop = " << loop+1 << "; Explained Variance (%) = " << (f/st)*100 << "; iter = " << it+1 << "; fdif = " << fdif << std::endl;
      if(loop==0){
        Vbest = V;
        Abest = A;
        Ybest = Xs*Abest;
        fbest = f;
        loopbest = 1;
        itbest = it+1;
        fdifbest = fdif;
      }
      if(f > fbest){
        Vbest = V;
        Abest = A;
        Ybest = Xs*A;
        fbest = f;
        loopbest = loop+1;
        itbest = it+1;
        fdifbest = fdif;
      }
    }
    // sort clusters of variables in descending order of variance
    varY = arma::var(Ybest, 1).t();
    ic = arma::sort_index(varY, "descend");
    varY = varY(ic);
    Abest = Abest.cols(ic);
    Vbest = Vbest.cols(ic);
    Ybest = Ybest.cols(ic);
    if(verbose==1)
      Rcpp::Rcout << "DPCA (Final): Percentage of Explained Variance (%) = " << (fbest/st)*100 << "; loop = " << loopbest << "; iter = " << itbest << "; fdif = " << fdifbest << std::endl;
    if(print==1){
      Rcpp::Rcout << "\n>> Variance explained by the DPCA (% BSS / TSS)= "<<(fbest/st)*100 << "\n" << std::endl;
      Rcpp::Rcout << ">> Loading Matrix (Manifest Variables x Latent variables) \n" << std::endl;
      Rcpp::CharacterVector namesc(Q);
      Rcpp::CharacterVector namesq = Nomi[1];
      for(int j = 0; j<Q; j++){
        namesc(j) = "PC " + Rcpp::toString(j+1);
      }
      Rcpp::NumericMatrix matrice = wrap(Abest);
      Rcpp::NumericMatrix Out(Q,6);
      Rcpp::colnames(matrice) = namesc;
      Rcpp::rownames(matrice) = namesq;
      Rf_eval(Rf_lang2(Rf_install("print"), matrice), R_GlobalEnv);
      
        arma::vec e2k(Q);
        arma::vec cro(Q);
        for(g = 0; g<Q; g++){
          ibCg = Vbest.col(g);
          if(arma::sum(ibCg)>1){
            indices = arma::find(ibCg==1);
            JCg = JJ.elem(indices);
            cro(g) = CronbachAlpha(Xs.cols(arma::conv_to<arma::uvec>::from(JCg)));
            Sg = S.submat(arma::conv_to<arma::uvec>::from(JCg), arma::conv_to<arma::uvec>::from(JCg));
            if(n > 20000)
              arma::eig_sym(l, U, Sg, "dc");
            else
              arma::eig_sym(l, U, Sg);
            dL = arma::sort(l, "descend");
            idL = arma::sort_index(l, "descend");
            L = arma::diagmat(dL);
            U = U.cols(idL);
            e2k(g)=L(1,1);
          }
          else{
            e2k(g) = 0;
            cro(g) = 1;
          }
        }
        Rcpp::Rcout << "\n>> Summary of the latent factors:" << std::endl;
        for(int q = 0; q<Q; q++){
          varcum = varcum + varY(q);
          Rcpp::NumericVector vals = {varY(q), (varY(q)/J)*100, varcum, (varcum/J)*100, e2k(q), cro(q)};
          Out(q,_) = vals;
        }
        Rcpp::rownames(Out) = namesc;
        
        Rcpp::CharacterVector fields= {"Explained Variance", "Expl. Var. (%)", "Cumulated Var.", "Cum. Var (%)", "Var. 2nd component", "Cronbach's Alpha"}; 
        Rcpp::colnames(Out) = fields;
        Function print("print");
        print(Out, Named("quote") = false);
      }
    return Rcpp::List::create(Rcpp::Named("A") = Abest,
                            Rcpp::Named("V") = Vbest,
                            Rcpp::Named("betweenss") = fbest,
                            Rcpp::Named("totss") = st,
                            Rcpp::Named("size") = arma::sum(Vbest,0),
                            Rcpp::Named("loop") = loopbest,
                            Rcpp::Named("it") = itbest);
  }
  else
    stop("Re-run with proper values for the arguments.");
}
