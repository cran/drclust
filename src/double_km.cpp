#define ARMA_WARN_LEVEL 0
#include "rc_aux.h"
#include "prep.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

//' @name doublekm
//' @title Double k-means Clustering
//' @description
//' Performs simultaneous \emph{k}-means partitioning on units and variables (rows and columns of the data matrix). 
//' 
//' @usage doublekm(Xs, K, Q, Rndstart, verbose, maxiter, tol, prep, print)
//' 
//' @param Xs Units x variables numeric data matrix.
//' @param K Number of clusters for the units.
//' @param Q Number of clusters for the variables.
//' @param Rndstart Number of runs to be performed (Defaults is 20).
//' @param verbose Outputs basic summary statistics for each run (1 = enabled; 0 = disabled, default option).
//' @param maxiter Maximum number of iterations allowed (if convergence is not yet reached. Default is 100).
//' @param tol Tolerance threshold. It is the maximum difference between the values of the objective function of two consecutive iterations such that convergence is assumed (default is 1e-6).
//' @param prep Pre-processing of the data. 1 performs the z-score transform (default choice); 2 performs the min-max transform; 0 leaves the data un-pre-processed.
//' @param print Prints summary statistics of the results (1 = enabled; 0 = disabled, default option).
//' 
//' @return returns a list of estimates and some descriptive quantities of the final results.
//' @return \item{U}{Units x clusters membership matrix (binary and row-stochastic). Each row is a dummy variable indicating to which unit-cluster each unit has been assigned.}
//' @return \item{V}{Variables x clusters membership matrix (binary and row-stochastic). Each row is a dummy variable indicating to which variable-cluster each variable has been assigned.}
//' @return \item{centers}{K x Q matrix of centers containing the row means expressed in terms of column means.}
//' @return \item{totss}{The total sum of squares (scalar).}
//' @return \item{withinss}{Vector of within-row-cluster sum of squares, one component per cluster.}
//' @return \item{columnwise_withinss}{Vector of within-column-cluster sum of squares, one component per cluster.}
//' @return \item{betweenss}{Amount of deviance captured by the model (scalar).}
//' @return \item{K-size}{Number of units assigned to each row-cluster (vector).}
//' @return \item{Q-size}{Number of variables assigned to each column-cluster (vector).}
//' @return \item{pseudoF}{Calinski-Harabasz index of the resulting (row-) partition (scalar).}
//' @return \item{loop}{The index of the (best) run from which the results have been chosen.}
//' @return \item{it}{the number of iterations performed during the (best) run.}
//' 
//' @author Ionel Prunila, Maurizio Vichi
//' 
//' @references Vichi M. (2001) "Double k-means Clustering for Simultaneous Classification of Objects and Variables" <doi:10.1007/978-3-642-59471-7_6>
//' 
//' @examples
//' # Iris data 
//' # Loading the numeric variables of iris data
//' iris <- as.matrix(iris[,-5]) 
//' 
//' # double k-means with 3 unit-clusters and 2 variable-clusters
//' out <- doublekm(iris, K = 3, Q = 2)
//' 
//' @export

//[[Rcpp::export]]

List doublekm(arma::mat Xs, int K, int Q, int Rndstart = 20, int verbose = 0, int maxiter = 100, double tol = 1e-6, int prep = 1, int print = 0){
  int n = Xs.n_rows, J = Xs.n_cols;
  arma::mat Xs2,  U, V, Ym, B, Ymv, Ymu, Ym2, cYm, cYm2, M, S2x, S1x;
  arma::mat Ymean2;
  int it, loop;
  double f0, fdif, f=0;
  arma::uvec iicc, ic;
  arma::rowvec sv, su, wd(K), cwd(Q);
  // best estimates
  arma::mat Vbest, Ubest, Abest, Ymbest, Ybest;
  double fbest=0, fdifbest=0, st;
  int loopbest, itbest;
  arma::rowvec wdbest(K), cwdbest(Q);

  bool args_ok = checkArgs(Q, Rndstart, verbose, maxiter, tol, J);
  bool k_ok = checkK(K, n);
  Xs = preproc(Xs, prep);
  bool prep_ok = checkPrep(prep);
  bool stats_ok = checkStats(print);
  
  Xs2 = arma::pow(Xs,2);
  st = sum(sum(Xs2));

  S2x = arma::repmat(arma::sum(arma::pow(Xs,2),1),1,K); // assign function
  S1x = arma::repmat(arma::sum(arma::pow(Xs.t(),2),1),1,Q);

  if(args_ok==1 && k_ok==1 && prep_ok == 1 && stats_ok == 1){
  for(loop = 0; loop < Rndstart; loop++){
    // initialization
    V = randPU(J,Q);
    U = randPU(n,K);
    su = arma::sum(U,0);
    sv = arma::sum(V,0);
    Ym = arma::diagmat(1.0/su)*U.t()*Xs*V*arma::diagmat(1.0/sv);
    B = U*Ym*arma::trans(V);
    f0 = arma::trace(arma::trans(B)*B);
    it = 0;
    fdif = 2*tol;
    //iteration phase
    while((fdif > tol) || (it >= maxiter)){
      it++;
      //Given Xm and V, update U
      Ymv = Ym*arma::trans(V);
      U = assign_ssed(S2x, Ymv, Xs, K, n);
      su = arma::sum(U,0);
      while(arma::sum(su==0)>0){ // Filling empty clusters
        Ym = arma::diagmat(1.0/su)*U.t()*Xs*V*arma::diagmat(1.0/sv);
        Ym2 = arma::diagmat(1.0/su)*U.t()*Xs2*V*arma::diagmat(1.0/sv);
        wd = arma::trans(arma::sum((arma::diagmat(su)*(Ym2 - arma::pow(Ym,2)) * V.t() ),1));
        wd.elem(arma::find(su==0 || su==1)).fill(0);
        U = split_maxwd(su, wd, U, Xs, K);
        su = arma::sum(U,0);
      }
      Ym = arma::diagmat(1.0/su)*U.t()*Xs*V*arma::diagmat(1.0/sv);

      //Given Xm and U, update V
      Ymu = U*Ym;
      M = Ym.t()*U.t();
      V = assign_ssed(S1x, M, Xs.t(), Q, J);
      sv = arma::sum(V,0);
      while(arma::sum(sv==0)>0){ // If verified, solve the empty cluster issue
        cYm = arma::diagmat(1/sv)*V.t()*Xs.t()*U*arma::diagmat(1/su);
        cYm2 = arma::diagmat(1/sv)*V.t()*Xs2.t()*U*arma::diagmat(1/su);
        cwd = arma::trans(arma::sum((arma::diagmat(sv)*(cYm2 - arma::pow(cYm,2))*U.t()),1));
        cwd.elem(arma::find(sv==0 || sv == 1)).fill(0);
        V = split_maxwd(sv, cwd, V, Xs.t(), Q);
        sv = arma::sum(V,0);
      }
      //within deviances
      Ym = arma::diagmat(1.0/su)*U.t()*Xs*V*arma::diagmat(1.0/sv);
      Ym2 = arma::diagmat(1.0/su)*U.t()*Xs2*V*arma::diagmat(1.0/sv);
      wd = arma::trans(arma::sum((arma::diagmat(su)*(Ym2 - arma::pow(Ym,2)) * V.t() ),1));

      cYm = arma::diagmat(1/sv)*V.t()*Xs.t()*U*arma::diagmat(1/su);
      cYm2 = arma::diagmat(1/sv)*V.t()*Xs2.t()*U*arma::diagmat(1/su);
      cwd = arma::trans(arma::sum((arma::diagmat(sv)*(cYm2 - arma::pow(cYm,2))*U.t()),1));
      // update f
      B = U*Ym*arma::trans(V);
      f=arma::trace(arma::trans(B)*B);
      fdif = f - f0;
      if(fdif>tol){
        f0=f;
      }
      else{
        break;
      }
    }
    if(verbose==1)
      Rcpp::Rcout <<"DKM: Loop: = " << loop + 1 << "; Explained variance (%) = " << (f/st)*100 << "; iter = " << it +1 << "; fdif: " << fdif << std::endl;
    //update estimates
    if(loop==0){
      Ubest=U;
      Vbest = V;
      Ymbest = Ym;
      fbest=f;
      fdifbest = fdif;
      loopbest = 1;
      itbest = it+1;
      wdbest = wd;
      cwdbest = cwd;
    }
    if(f>fbest){
      Ubest=U;
      Vbest = V;
      Ymbest = Ym;
      fbest = f;
      fdifbest = fdif;
      loopbest = loop+1;
      itbest = it+1;
      wdbest = wd;
      cwdbest = cwd;
    }
  }
  //sort clusters of variables per descending order of cardinality
  ic = arma::sort_index(arma::diagvec(Vbest.t()*Vbest), "descend");
  Vbest = Vbest.cols(ic);
  cwdbest = arma::trans(cwdbest(ic));

  //sort clusters of objects in descending order of cardinality
  ic = arma::sort_index(arma::diagvec(Ubest.t()*Ubest), "descend");
  Ubest = Ubest.cols(ic);
  wdbest = arma::trans(wdbest(ic));
  double pseudoF = (f/(K-1))/((st-f)/(n-K));
  if(verbose==1)
    Rcpp::Rcout <<"DKM (Final) Explained variance (%) = " << (fbest/st)*100 << "; loop = " << loopbest  << "; iter = " << itbest << "; fidf = " << fdifbest << std::endl;
  if(print==1){
    Rcpp::Rcout << "\n>> Variance Explained by the DKM (% BSS / TSS):  " << (fbest/st)*100 << std::endl;
    Rcpp::Rcout << "\n>> Centroid Matrix (Unit-centroids x Variable-centroids):\n" << std::endl;
    Rcpp::NumericMatrix centroids = wrap(Ymbest);
    Rcpp::CharacterVector namesr(K);
    Rcpp::CharacterVector namesc(Q);
    
    Rcpp::NumericVector rowdev = wrap(wdbest);
    Rcpp::NumericVector coldev = wrap(cwdbest);
    Rcpp::NumericVector rowsize = trunc(wrap(sum(Ubest,0)));
    Rcpp::NumericVector colsize = trunc(wrap(sum(Vbest,0)));
    Rcpp::CharacterMatrix Outu(2,K);
    Rcpp::CharacterMatrix Outv(2,Q);
    
    for(int i = 0; i<K; i++){
      namesr(i) = "U-Clust " + Rcpp::toString(i+1);
      Outu(0,i) = rowsize(i);
      Outu(1,i) = rowdev(i);
    }
    for(int j = 0; j<Q; j++){
      namesc(j) = "V-Clust " + Rcpp::toString(j+1);
      Outv(0,j) = colsize(j);
      Outv(1,j) = coldev(j);
    }
    Rcpp::rownames(centroids) = namesr;
    Rcpp::colnames(centroids) = namesc;
    
    Rf_eval(Rf_lang2(Rf_install("print"), centroids), R_GlobalEnv);
    
    Rcpp::Rcout << "\n>> Unit-clusters: \n" << std::endl;
    Rcpp::colnames(Outu) = namesr;
    Rcpp::StringVector fields = {"Size", "Deviance"};
    Rcpp::rownames(Outu) = fields;
    Function print("print");
    print(Outu, Named("quote") = false);
    
    Rcpp::Rcout << "\n>> Variable-clusters: \n " << std::endl;
    Rcpp::colnames(Outv) = namesc;
    Rcpp::rownames(Outv) = fields;
    print(Outv, Named("quote") = false);
    
    Rcpp::Rcout << "\n>> pseudoF Statistic (Calinski-Harabasz): " << pseudoF << std::endl;
  }
  return Rcpp::List::create(Rcpp::Named("U") = Ubest,
                              Rcpp::Named("V") = Vbest,
                               Rcpp::Named("centers") = Ymbest,
                               Rcpp::Named("withinss") = wdbest,
                               Rcpp::Named("columnwise_withinss") = cwdbest,
                               Rcpp::Named("betweenss") = fbest,
                               Rcpp::Named("totss") = st,
                               Rcpp::Named("K-size") = sum(Ubest,0),
                               Rcpp::Named("Q-size") = sum(Vbest,0),
                               Rcpp::Named("pseudoF") = pseudoF,
                               Rcpp::Named("loop") = loopbest,
                               Rcpp::Named("it") = itbest);
  }
  else
    stop("Re-run with proper values for the arguments.");
}
