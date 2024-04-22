#define ARMA_WARN_LEVEL 0
#include "rc_aux.h"
#include "dispca.h"
#include "prep.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]


//' @name dpcakm
//' @title Clustering with Disjoint Principal Components Analysis
//' @description
//' Performs simultaneously k-means partitioning on units and disjoint PCA on the variables, computing each principal component from a different subset of variables. The result is a simplified, easier to interpret loading matrix A, 
//' the principal components and the clustering. The reduced subspace is identified by the centroids.
//' 
//' 
//' @usage dpcakm(X, K, Q, Rndstart, verbose, maxiter, tol, constr, print, prep)
//' 
//' @param X Units x variables numeric data matrix.
//' @param K Number of clusters for the units.
//' @param Q Number of principal components.
//' @param Rndstart Number of runs to be performed (Defaults is 20).
//' @param verbose Outputs basic summary statistics for each run (1 = enabled; 0 = disabled, default option).
//' @param maxiter Maximum number of iterations allowed (if convergence is not yet reached. Default is 100).
//' @param tol Tolerance threshold (maximum difference between the values of the objective function of two consecutive iterations such that convergence is assumed. Default is 1e-6).
//' @param constr is a vector of length J = nr. of variables, pre-specifying to which cluster some of the variables must be assigned. Each component of the vector can assume integer values from 1 o Q = nr. of variable-cluster / principal components (See examples for more details), or 0 if no constraint on the variable is imposed (i.e., it will be assigned based on the plain algorithm).
//' @param print Prints summary statistics of the results (1 = enabled; 0 = disabled, default option).
//' @param prep Pre-processing of the data. 1 performs the z-score transform (default choice); 2 performs the min-max transform; 0 leaves the data un-pre-processed.
//'  
//' 
//' @return returns a list of estimates and some descriptive quantities of the final results.
//' @return \item{V}{Variables x factors membership matrix (binary and row-stochastic). Each row is a dummy variable indicating to which cluster each variable has been assigned.}
//' @return \item{U}{Units x clusters membership matrix (binary and row-stochastic). Each row is a dummy variable indicating to which cluster each unit has been assigned.}
//' @return \item{A}{Variables x components loading matrix.}
//' @return \item{centers}{K x Q matrix of centers containing the row means expressed in the reduced space of Q principal components.}
//' @return \item{totss}{The total sum of squares (scalar).}
//' @return \item{withinss}{Vector of within-cluster sum of squares, one component per cluster.}
//' @return \item{betweenss}{Amount of deviance captured by the model (scalar).}
//' @return \item{K-size}{Number of units assigned to each row-cluster (vector).}
//' @return \item{Q-size}{Number of variables assigned to each column-cluster (vector).}
//' @return \item{pseudoF}{Calinski-Harabasz index of the resulting partition (scalar).}
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
//' out <- dpcakm(iris, K = 3, Q = 2, Rndstart = 5)
//' 
//' # Constraint: the first two variables must contribute to the same factor.
//' outc <- dpcakm(iris, K = 3, Q = 2, Rndstart = 5,constr = c(1,1,0,0))
//' @export
 


//[[Rcpp::export]]

List dpcakm(Rcpp::NumericMatrix X, int K, int Q, int Rndstart = 20, int verbose = 0, int maxiter = 100, double tol = 1e-6, arma::vec constr = 00, int print = 0, int prep = 1){
  Rcpp::List Nomi = X.attr("dimnames");
  arma::mat Xs = as<arma::mat>(X);
  
  /// Declaration of variables
  int n = Xs.n_rows, J = Xs.n_cols;
  arma::mat Xs2, S, U, V, Ym, Xmean, XX, S2x;
  arma::mat L, A, A0, Y, Ymean, Xmean2, Ymean2, Sg, Av, VC = arma::eye(Q,Q), S0g, AAv;
  int it, loop, j, g, a, aa, posmax, flg;
  double f0, fdif, f=0, fmax, st, pseudoF;
  arma::vec un, l, ll, dL, varY, JJ, ibCg, JCg, ibCpm, JCpm, y, zJ = arma::zeros(J,1);
  arma::uvec idL, iicc, indices, ic;
  arma::rowvec su, wd(K);
  
  // best estimates
  arma::mat Vbest, Ubest, Abest, Ymbest, Ybest, Xmeanbest;
  double fbest=0, fdifbest=0;
  int loopbest, itbest;
  arma::rowvec wdbest(K);
  
  
  Xs = preproc(Xs, prep);
  bool prep_ok = checkPrep(prep);
  bool stats_ok = checkStats(print);
  Xs = Xs*(n/(n-1));
  Xs2 = pow(Xs,2);
  st = sum(sum(Xs2));
  
  S = cov(Xs,1);
  S2x = repmat(sum(pow(Xs,2),1),1,K); // assign function
  bool args_ok = checkArgs(Q, Rndstart, verbose, maxiter, tol, J);
  bool k_ok = checkK(K, n);
  bool constr_ok = checkConstr(constr, J, Q);
  if(constr.size()==1)
    constr = arma::zeros(J,1);
  JJ = arma::regspace(0,1,J-1);
  if(args_ok == 1 && k_ok == 1 && constr_ok ==1 && stats_ok ==1 && prep_ok == 1){
    for(loop = 0; loop < Rndstart; loop++){
      // Initialization
      flg = 1;
      //////////////it = 0;
      while(flg>0){ //constraint -> CFA
        V = randPU(J,Q);
        for(j=0;j<J;j++){
          if(constr(j)>0){
            V.row(j) = VC.row(constr(j)-1);
          }
        }
        flg = arma::sum(arma::find(arma::sum(V,0)==0));
      }
      U = randPU(n,K);
      su = arma::sum(U,0);
      A = arma::zeros(J,Q);
      //update centroid matrix Xmean
      Xmean = arma::diagmat(1.0/su)*U.t()*Xs;
      
      //Rcpp::Rcout<<"fino a qui tutto bene" << std::endl;
       
      for(g = 0; g < Q; g++){
        ibCg = V.col(g);
        JCg = JJ.elem(arma::find(ibCg==1));
        //Sg = S.submat(arma::conv_to<arma::uvec>::from(JCg), arma::conv_to<arma::uvec>::from(JCg));
        S = Xs.t()*U*arma::diagmat(1.0/su)*U.t()*Xs;
        if(arma::sum(ibCg) > 1){
          Sg = S.submat(arma::conv_to<arma::uvec>::from(JCg), arma::conv_to<arma::uvec>::from(JCg));
          if(n>10000)
            arma::eig_sym(l, Av, Sg, "dc");
          else
            arma::eig_sym(l, Av, Sg);
          a = arma::index_max(abs(l));
          y = Xs.cols(arma::conv_to<arma::uvec>::from(JCg))*Av.col(a);
          A.submat(arma::conv_to<arma::uvec>::from(JCg),arma::uvec{static_cast<unsigned int>(g)}) = Av.col(a);
          
        }
        else{
          A.submat(arma::conv_to<arma::uvec>::from(JCg), arma::uvec{static_cast<unsigned int>(g)}) = arma::ones(JCg.n_elem);
        }
        
      }
      Ymean = Xmean*A;
      f0 = arma::trace(Ymean.t()*U.t()*U*Ymean);
      fmax = 0;
      //it = 0;
      fdif = 2*tol;
      //iteration phase
      //while(fdif > tol || it >= maxiter){
      for(it=0; it < maxiter; it++){
        //given Ymean and A, update U
        Y = Xs*A;
        S2x = arma::repmat(arma::sum(arma::pow(Y,2),1),1,K);
        U = assign_ssed(S2x, Ymean, Y, K, n);
        su = arma::sum(U,0);
        while(arma::sum(su==0)>0){
          Xmean = arma::diagmat(1.0/su)*U.t()*Xs;
          Xmean2 = arma::diagmat(1.0/su)*U.t()*Xs2;
          Ymean2 = arma::diagmat(1.0/su)*U.t()*arma::pow(Xs*A,2);
          wd = arma::trans(arma::sum((arma::diagmat(su)*(Ymean2 - arma::pow(Xmean*A,2))),1));
          wd.elem(arma::find(su==0 || su==1)).fill(0);
          U = split_maxwd(su, wd, U, Xs*A, K);
          su = arma::sum(U,0);
        }
        //Given U and A, compute Xmean (centroids)
        Xmean = diagmat(1.0/su)*U.t()*Xs;
        //Given U and Ymean, update V and A
        S = Xs.t()*U*arma::diagmat(1.0/su)*U.t()*Xs;
        for(j = 0; j < J; j++){
          posmax = JJ(arma::index_max(V.row(j)==1));
          if(constr(j)==0){
            for(g = 0; g < Q; g++){
              V.row(j) = VC.row(g);
              if(sum(V.col(posmax))>0){
              ibCg = V.col(g); // new class of V
              ibCpm = V.col(posmax); // old class of V
              JCg = JJ.elem(arma::find(ibCg==1));
              JCpm = JJ.elem(arma::find(ibCpm==1));
              
              if(arma::sum(ibCg)>1){
                Sg = S.submat(arma::conv_to<arma::uvec>::from(JCg),arma::conv_to<arma::uvec>::from(JCg));
                if(n>10000)
                  arma::eig_sym(l, Av, Sg, "dc");
                else
                  arma::eig_sym(l, Av, Sg);
                a = arma::index_max(abs(l));
                y = Xs.cols(arma::conv_to<arma::uvec>::from(JCg))*Av.col(a);
                
                if(arma::sum(Av.col(a))<0){
                  Av.col(a) = -Av.col(a);
                }
                A.col(g) = zJ;
                
                A.submat(arma::conv_to<arma::uvec>::from(JCg), arma::uvec{static_cast<unsigned int>(g)}) = Av.col(a);
                
              }
              else{
                A.col(g) = zJ;
                A.submat(arma::conv_to<arma::uvec>::from(JCg), arma::uvec{static_cast<unsigned int>(g)}) = arma::ones(JCg.n_elem);
              }
              
              if(sum(ibCpm)>1){
                Sg=S.submat(arma::conv_to<arma::uvec>::from(JCpm),arma::conv_to<arma::uvec>::from(JCpm));
                
                if(n>10000)
                  arma::eig_sym(ll, AAv, Sg, "dc");
                else
                  arma::eig_sym(ll, AAv, Sg);
                aa = arma::index_max(abs(ll));
                y = Xs.cols(arma::conv_to<arma::uvec>::from(JCpm))*AAv.col(aa);
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
                
              //A = dpca_updateA(posmax, g, ibCg, ibCpm, S, A, JJ, J);
              Ymean = Xmean*A;
              f = arma::trace(Ymean.t()*(U.t()*U)*Ymean);
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
        //within sum of squares
        Ymean = Xmean*A;
        Ymean2 = arma::diagmat(1.0/su)*U.t()*arma::pow(Xs*A,2);
        Xmean2 = arma::diagmat(1.0/su)*U.t()*Xs2;
        wd = arma::trans(arma::sum((arma::diagmat(su)*(Ymean2 - arma::pow(Xmean*A,2))),1));
        // f
        f = arma::trace(Ymean.t()*(U.t()*U)*Ymean);
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
      Ymean = Xmean*A;
      f = arma::trace(Ymean.t()*(U.t()*U)*Ymean);
      fdif = f - f0;
      if(verbose==1)
        Rcpp::Rcout << "CDPCA: Loop = " << loop+1 << "; Explained Variance (%) = " << (f/st)*100 << "; iter = " << it+1 << "; fdif = " << fdif << std::endl;
      if(loop==0){
        Vbest = V;
        Ubest = U;
        Abest = A;
        Ybest = Xs*Abest;
        fbest = f;
        Ymbest = Ymean;
        loopbest = 1;
        itbest = it+1;
        fdifbest = fdif;
        wdbest = wd;
      }
      if(f > fbest){
        Vbest = V;
        Ubest = U;
        fbest = f;
        Abest = A;
        Ybest = Xs*Abest;
        Ymbest = Ymean;
        loopbest = loop+1;
        itbest = it+1;
        fdifbest = fdif;
        wdbest = wd;
      }
    }
    Ybest = Xs*Abest;
    // sort clusters of variables in descending order of variance
    varY = arma::trans(arma::var(Ybest, 1));
    ic = arma::sort_index(varY, "descend");
    Abest = Abest.cols(ic);
    Vbest = Vbest.cols(ic);
    Ybest = Ybest.cols(ic);
    // sort clusters of objects in descending order of variance
    ic = arma::sort_index(arma::diagvec(Ubest.t()*Ubest), "descend");
    Ubest = Ubest.cols(ic);
    wdbest = arma::trans(wdbest(ic));
    pseudoF = (f/(K-1))/((st-f)/(n-K));
    if(verbose==1)
      Rcpp::Rcout << "CDPCA (Final): Percentage of Explained Variance (%) = " << (fbest/st)*100 << "; loop = " << loopbest << "; iter = " << itbest << "; fdif = " << fdifbest << std::endl;
    if(print==1){
      Rcpp::Rcout << "\n>> Variance Explained by the DPCAKM (% BSS / TSS): " << (fbest/st)*100 << std::endl;
      
      Rcpp::Rcout << "\n>> Matrix of Centroids (Unit-centroids x Principal Components):\n" << std::endl;
      
      Rcpp::NumericMatrix matrice = wrap(Ymbest);
      Rcpp::CharacterVector namesr(K);
      Rcpp::CharacterVector namesc(Q);
      Rcpp::CharacterVector namesq = Nomi[1];
      Rcpp::NumericVector rowdev = wrap(wdbest);
      
      Rcpp::NumericVector rowsize = trunc(wrap(sum(Ubest,0)));
      Rcpp::NumericVector colsize = trunc(wrap(sum(Vbest,0)));
      Rcpp::CharacterMatrix Outu(2,K);
      Rcpp::CharacterMatrix Outq(Q,4);

      for(int i = 0; i<K; i++){
        namesr(i) = "Clust " + Rcpp::toString(i+1);
        Outu(0,i) = rowsize(i);
        Outu(1,i) = rowdev(i);
      }
      for(int j = 0; j<Q; j++){
        namesc(j) = "PC " + Rcpp::toString(j+1);
      }
      Rcpp::rownames(matrice) = namesr;
      Rcpp::colnames(matrice) = namesc;
      Rf_eval(Rf_lang2(Rf_install("print"), matrice), R_GlobalEnv);
      
      Rcpp::colnames(Outu) = namesr;
      Rcpp::StringVector fields = {"Size", "Deviance"};
      Rcpp::rownames(Outu) = fields;
      Rcpp::Rcout << "\n>> Unit-clusters: " << std::endl;
      Function print("print");
      print(Outu, Named("quote") = false);
      Rcpp::Rcout << "\n>> Loading Matrix (Manifest Variables x Latent Variables):\n " << std::endl;
      matrice = wrap(Abest);
      Rcpp::colnames(matrice) = namesc;
      Rcpp::rownames(matrice) = namesq;
      Rf_eval(Rf_lang2(Rf_install("print"), matrice), R_GlobalEnv);
      
      double varcum = 0;
      Rcpp::Rcout << "\n>> Summary of the latent factors:" << std::endl;
      for(int q = 0; q<Q; q++){
        varcum = varcum + varY(q);
        Rcpp::NumericVector vals = {varY(q), (varY(q)/J)*100, varcum, (varcum/J)*100};
        Outq(q,_) = vals;
      }
      Rcpp::rownames(Outq) = namesc;
      
      fields= {"Explained Variance", "Expl. Var. (%)", "Cumulated Var.", "Cum. Var (%)"}; 
      Rcpp::colnames(Outq) = fields;
      print(Outq, Named("quote") = false);
      
      Rcpp::Rcout << "\n>> pseudoF Statistic (Calinski-Harabasz): " << pseudoF << std::endl;
    }
  return Rcpp::List::create(Rcpp::Named("U") = Ubest,
                            Rcpp::Named("A") = Abest,
                            Rcpp::Named("V") = Vbest,
                            Rcpp::Named("centers") = Ymbest,
                            Rcpp::Named("withinss") = wdbest,
                            Rcpp::Named("betweenss") = fbest,
                            Rcpp::Named("totss") = st,
                            Rcpp::Named("K-size") = arma::sum(Ubest,0),
                            Rcpp::Named("Q-size") = arma::sum(Vbest,0),
                            Rcpp::Named("pseudoF") = pseudoF,
                            Rcpp::Named("loop") = loopbest,
                            Rcpp::Named("it") = itbest);
  }
  else
    stop("Re-run with proper values for the arguments.");
}
