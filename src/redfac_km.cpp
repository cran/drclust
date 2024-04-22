#define ARMA_WARN_LEVEL 0
#include "rc_aux.h"
#include "prep.h"
#include <RcppArmadillo.h>


using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

arma::mat varimax(arma::mat A){

  double conv = 0.000001;
  int m = A.n_rows;
  int r = A.n_cols;
  int i, j, sign;
  int iter = 0;
  arma::mat T = arma::eye(r, r);
  arma::mat B = A;
  arma::vec x, y, xx, yy, u, v, w, vv, ww;
  double a, b, c, sin, cos, vvv=0;
  arma::vec mones = arma::ones(m,1);

  double f = arma::accu(arma::pow((arma::pow(A,2) - mones*arma::sum(arma::pow(A,2),0)/m),2));
  double fold = f-2*conv*f;
  if(f==0){
    fold = -conv;
  }
  while(f-fold > f*conv){
    fold = f;
    iter++;
    for(i = 0; i<r; i++){
      for(j = i+1; j<r; j++){
        x = B.col(i);
        y = B.col(j);
        xx = T.col(i);
        yy = T.col(j);
        u = arma::square(x)-arma::square(y);
        v = 2*x%y;
        u = u - mones*arma::sum(u)/m;
        v = v - mones*arma::sum(v)/m;
        a = 2*arma::sum(u%v);
        b = arma::sum(arma::square(u))-arma::sum(arma::square(v));
        c = std::sqrt((a*a) + (b*b));
        if(a>=0){
          sign = 1;
        }
        //if(a < 0){ 
        else{
          sign = -1;
        }
        if(c < 0.00000000001){
          cos = 1;
          sin = 0;
        }
        else{
          vvv = -sign * std::sqrt((b+c)/(2*c));
          sin = std::sqrt(0.5 - 0.5*vvv);
          cos = std::sqrt(0.5 + 0.5*vvv);
        }
        v = cos*x - sin*y;
        w = cos*y + sin*x;
        vv = cos*xx - sin*yy;
        ww = cos*yy + sin*xx;
        if(vvv>=0){
          B.col(i) = v;
          B.col(j) = w;
          T.col(i) = vv;
          T.col(j) = ww;
        }
        else{
          B.col(j) = v;
          B.col(i) = w;
          T.col(j) = vv;
          T.col(i) = ww;
        }
      }
    }
    f = arma::accu(arma::pow((arma::pow(B,2) - mones*arma::sum(arma::pow(B,2),0)/m),2));
  }
  return B;
}


//' @name redkm
//' @title k-means on a reduced subspace
//' @description
//' Performs simultaneously k-means partitioning on units and principal component analysis on the variables. 
//' 
//' @usage redkm(X, K, Q, Rndstart, verbose, maxiter, tol, rot, prep, print)
//' 
//' @param X Units x variables numeric data matrix.
//' @param K Number of clusters for the units.
//' @param Q Number of principal components w.r.t. variables.
//' @param Rndstart Number of runs to be performed (Defaults is 20).
//' @param verbose Outputs basic summary statistics for each run (1 = enabled; 0 = disabled, default option).
//' @param maxiter Maximum number of iterations allowed (if convergence is not yet reached. Default is 100).
//' @param tol Tolerance threshold (maximum difference between the values of the objective function of two consecutive iterations such that convergence is assumed. Default is 1e-6).
//' @param rot performs varimax rotation of axes obtained via PCA. (=1 enabled; =0 disabled, default option)
//' @param prep Pre-processing of the data. 1 performs the z-score transform (default choice); 2 performs the min-max transform; 0 leaves the data un-pre-processed.
//' @param print Tolerancestats summary statistics of the performed method (1 = enabled; 0 = disabled, default option).
//' 
//' 
//' 
//' @return returns a list of estimates and some descriptive quantities of the final results.
//' @return \item{U}{Units x clusters membership matrix (binary and row-stochastic). Each row is a dummy variable indicating to which cluster each unit has been assigned.}
//' @return \item{A}{Variables x components loading matrix (orthonormal).}
//' @return \item{centers}{K x Q matrix of centers containing the row means expressed in the reduced space of Q principal components.}
//' @return \item{totss}{The total sum of squares (scalar).}
//' @return \item{withinss}{Vector of within-cluster sum of squares, one component per cluster.}
//' @return \item{betweenss}{Amount of deviance captured by the model (scalar).}
//' @return \item{size}{Number of units assigned to each cluster (vector).}
//' @return \item{pseudoF}{Calinski-Harabasz index of the resulting partition (scalar).}
//' @return \item{loop}{The index of the (best) run from which the results have been chosen.}
//' @return \item{it}{the number of iterations performed during the (best) run.}
//' 
//' @author Ionel Prunila, Maurizio Vichi
//' 
//' @references 
//' de Soete G., Carroll J. (1994) "K-means clustering in a low-dimensional Euclidean space" <doi:10.1007/978-3-642-51175-2_24>
//' 
//' Kaiser H.F. (1958) "The varimax criterion for analytic rotation in factor analysis" <doi:10.1007/BF02289233>
//' 
//' @examples
//' # Iris data 
//' # Loading the numeric variables of iris data
//' iris <- as.matrix(iris[,-5]) 
//' 
//' # reduced k-means with 3 unit-clusters and 2 components for the variables
//' out <- redkm(iris, K = 3, Q = 2, Rndstart = 15, verbose = 0, maxiter = 100, tol = 1e-7, rot = 1)
//' 
//' @export
 


//[[Rcpp::export]]

List redkm(Rcpp::NumericMatrix X, int K, int Q, int Rndstart = 20, int verbose = 0, int maxiter = 100, double tol = 1e-6, int rot = 0, int prep = 1, int print = 0){
  Rcpp::List Nomi = X.attr("dimnames");
  arma::mat Xs = as<arma::mat>(X);
  /// Declaration of variables
  int n = Xs.n_rows, J = Xs.n_cols;
  arma::mat  Xs2, S, U, V, Ym, B, Ym2, Xmean, XX, S2x;
  arma::mat L, A, A0, Y, Ymean, RotA, Xmean2, Ymean2, Ymeanbest;
  int it, loop;
  double f0, fdif, f=0, st;
  arma::vec l, dL, varY;
  arma::uvec idL, iicc, ic;
  arma::rowvec su, wd(K);
  
  // best estimates
  arma::mat Vbest, Ubest, Abest, Ymbest, Ybest, Xmeanbest;
  double fbest=0, fdifbest=0;
  int loopbest, itbest;
  arma::rowvec wdbest(K);
  
  bool args_ok = checkArgs(Q, Rndstart, verbose, maxiter, tol, J);
  bool k_ok = checkK(K, n);
  bool rot_ok = 1;
  
  if(rot != 0 && rot != 1){
    Rcpp::Rcout << "Error: rot must be an integer = 0 or = 1" << std::endl;
    rot_ok = 0;
  }
  Xs = preproc(Xs, prep);
  bool prep_ok = checkPrep(prep);
  bool stats_ok = checkStats(print);
  
  Xs2 = arma::pow(Xs,2);
  st = arma::sum(arma::sum(Xs2));
  S = arma::cov(Xs,1);
  S2x = arma::repmat(arma::sum(arma::pow(Xs,2),1),1,K); // assign function
  if(args_ok ==1 && k_ok ==1 && rot_ok ==1 && prep_ok == 1 && stats_ok == 1){
    for(loop=0; loop<Rndstart; loop++){
      //Initialization
      U = randPU(n,K);
      su = arma::sum(U,0);
      Xmean = arma::diagmat(1.0/su)*U.t()*Xs;
      // Update A
      XX = Xs.t()*U*arma::diagmat(1.0/su)*U.t()*Xs;
      if(n>20000)
        arma::eig_sym(l, A, XX, "dc");
      else
        arma::eig_sym(l, A, XX);
      dL = arma::sort(l, "descend");
      idL = arma::sort_index(l, "descend");
      L = arma::diagmat(dL);
      A = A.cols(idL);
      A = A.cols(0,Q-1);
      // Update Ymean
      Ymean = Xmean*A;
      f0 = arma::trace(Ymean.t()*U.t()*U*Ymean);
      fdif=2*tol;
      it = 0;
      // iteration phase
      while(fdif>tol || it >= maxiter){ //
        it++;
        U = assign_ssed(S2x, Ymean*A.t(), Xs, K, n);
        su = arma::sum(U,0);
        while(arma::sum(su==0)>0){ // If verified, solve the empty cluster issue
          Xmean = arma::diagmat(1/su)*U.t()*Xs;
          Xmean2 = arma::diagmat(1/su)*U.t()*Xs2;
          wd = arma::trans(arma::sum((arma::diagmat(su)*(Xmean2 - arma::pow(Xmean*A*A.t(),2))),1));
          wd.elem(arma::find(su==0 || su==1)).fill(0);
          U = split_maxwd(su, wd, U, Xs, K);
          su = arma::sum(U,0);
        }
        Xmean = arma::diagmat(1.0/su)*U.t()*Xs;
        // Given Y and Xmean, update A
        XX = Xs.t()*U*arma::diagmat(1.0/su)*U.t()*Xs;
        if(n>20000)
          arma::eig_sym(l, A, XX, "dc");
        else
          arma::eig_sym(l, A, XX);
        dL = arma::sort(l, "descend");
        idL = arma::sort_index(l, "descend");
        L = arma::diagmat(dL);
        A = A.cols(idL);
        A = A.cols(0,Q-1);
        //Update Ymean and Y
        Ymean = Xmean*A;
        Y = Xs*A;
        // within sum of squares
        Xmean2 = arma::diagmat(1/su)*U.t()*Xs2;
        wd = arma::trans(arma::sum((arma::diagmat(su)*(Xmean2 - arma::pow(Xmean*A*A.t(),2))),1));
        f = arma::trace(Ymean.t()*U.t()*U*Ymean);
        fdif = f-f0;
        if(fdif>tol){
          f0=f;
          A0=A;
        }
        else{
          break;
        }
      }
      if(verbose==1)
        Rcpp::Rcout << "RKM: Loop = " << loop + 1 << "; Explained variance (%) = " << (f/st)*100 << "; iter = " << it+1 << "; fdif = " << fdif << std::endl;
      if(loop==0){
        Ubest=U;
        Abest=A;
        Ybest = Xs*Abest;
        Xmeanbest = Xmean;
        fbest = f;
        loopbest = 1;
        itbest = it+1;
        fdifbest = fdif;
        wdbest = wd;
      }
      if(f>fbest){
        Ubest=U;
        Abest = A;
        Ybest=Xs*Abest;
        Xmeanbest = Xmean;
        fbest=f;
        loopbest = loop+1;
        itbest = it+1;
        fdifbest = fdif;
        wdbest = wd;
      }
    }
    //sort components in descend order of variance
    varY = arma::trans(arma::var(Ybest,1));
    ic = arma::sort_index(varY, "descend");
    varY = varY(ic);
    Abest = Abest.cols(ic);
    if(rot != 0)
      Abest = varimax(Abest);
    Ybest = Ybest.cols(ic);
    Ymeanbest = Xmeanbest*Abest;
    iicc = arma::sort_index(arma::diagvec(Ubest.t()*Ubest), "descend");
    Ubest = Ubest.cols(iicc);
    wdbest = arma::trans(wdbest(iicc));
    double pseudoF = (f/(K-1))/((st-f)/(n-K));
    if(verbose==1)
      Rcpp::Rcout << "RKM (Final): Explained Variance (%) = " << (fbest/st)*100 << "; loop = "<< loopbest << "; iter = " << itbest << "; fdif = " << fdifbest << std::endl;
    if(print==1){
      Rcpp::Rcout << "\n>> Variance Explained by the RKM (% BSS / TSS): " << (fbest/st)*100 << std::endl;
      Rcpp::Rcout << "\n>> Matrix of Centroids (Unit-centroids x Principal Components):\n" << std::endl;
      
      Rcpp::NumericMatrix matrice = wrap(Ymeanbest);
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
      Rcpp::Rcout << "\n>> Summary of the latent factors:\n" << std::endl;
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
                            Rcpp::Named("centers") = Ymeanbest,
                            Rcpp::Named("withinss") = wdbest,
                            Rcpp::Named("betweenss") = fbest,
                            Rcpp::Named("totss") = st,
                            Rcpp::Named("size") = arma::sum(Ubest,0),
                            Rcpp::Named("pseudoF") = pseudoF,
                            Rcpp::Named("loop") = loopbest,
                            Rcpp::Named("it") = itbest);
  }
  else
    stop("Re-run with proper values for the arguments.");
}


//' @name factkm
//' @title Factorial k-means
//' @description
//' Performs simultaneously k-means partitioning on units and principal component analysis on the variables. 
//' Identifies the best partition in a Least-Squares sense in the best reduced space of the data. Both the data 
//' and the centroids are used to identify the best Least-Squares reduced subspace, where also their distances is measured.
//' 
//' 
//' @usage factkm(X, K, Q, Rndstart, verbose, maxiter, tol, rot, prep, print)
//' 
//' @param X Units x variables numeric data matrix.
//' @param K Number of clusters for the units.
//' @param Q Number of principal components w.r.t. variables.
//' @param Rndstart Number of runs to be performed (Defaults is 20).
//' @param verbose Outputs basic summary statistics for each run (1 = enabled; 0 = disabled, default option).
//' @param maxiter Maximum number of iterations allowed (if convergence is not yet reached. Default is 100).
//' @param tol Tolerance threshold (maximum difference in the values of the objective function of two consecutive iterations such that convergence is assumed. Default is 1e-6).
//' @param rot performs varimax rotation of axes obtained via PCA. (=1 enabled; =0 disabled, default option)
//' @param prep Pre-processing of the data. 1 performs the z-score transform (default choice); 2 performs the min-max transform; 0 leaves the data un-pre-processed.
//' @param print Prints summary statistics of the results (1 = enabled; 0 = disabled, default option).
//' 
//' 
//' @return returns a list of estimates and some descriptive quantities of the final results.
//' @return \item{U}{Units x clusters membership matrix (binary and row-stochastic). Each row is a dummy variable indicating to which cluster each unit has been assigned.}
//' @return \item{A}{Variables x components loading matrix (orthonormal).}
//' @return \item{centers}{K x Q matrix of centers containing the row means expressed in the reduced space of Q principal components.}
//' @return \item{totss}{The total sum of squares.}
//' @return \item{withinss}{Vector of within-cluster sum of squares, one component per cluster.}
//' @return \item{betweenss}{amount of deviance captured by the model.}
//' @return \item{size}{Number of units assigned to each cluster.}
//' @return \item{pseudoF}{Calinski-Harabasz index of the resulting partition.}
//' @return \item{loop}{The index of the (best) run from which the results have been chosen.}
//' @return \item{it}{the number of iterations performed during the (best) run.}
//' 
//' @author Ionel Prunila, Maurizio Vichi
//' 
//' @references 
//' Vichi M., Kiers H.A.L. (2001) "Factorial k-means analysis for two-way data" <doi:10.1016/S0167-9473(00)00064-5>
//' 
//' Kaiser H.F. (1958) "The varimax criterion for analytic rotation in factor analysis" <doi:10.1007/BF02289233>
//' 
//' @examples
//' # Iris data 
//' # Loading the numeric variables of iris data
//' iris <- as.matrix(iris[,-5]) 
//' 
//' # factorial k-means with 3 unit-clusters and 2 components for the variables
//' out <- factkm(iris, K = 3, Q = 2, Rndstart = 15, verbose = 0, maxiter = 100, tol = 1e-7, rot = 1)
//' 
//' @export


//[[Rcpp::export]]

List factkm(Rcpp::NumericMatrix X, int K, int Q, int Rndstart = 20, int verbose = 0, int maxiter = 100, double tol = 1e-6, int rot = 0, int prep = 1, int print = 0){
  Rcpp::List Nomi = X.attr("dimnames");
  arma::mat Xs = as<arma::mat>(X);
  /// Declaration of variables
  int n = Xs.n_rows, J = Xs.n_cols;
  arma::mat  Xs2, S, U, Ym, Xmean, XX, S2x;
  arma::mat L, A, A0, Y, Ymean, RotA, Xmean2, Ymean2;
  int it, loop;
  double f0, fdif, f=0, st;
  arma::vec un, l, dL, varY;
  arma::uvec idL, iicc, ic;
  arma::rowvec su, wd(K);
  
  // best estimates
  arma::mat Vbest, Ubest, Abest, Ymbest, Ybest, Xmeanbest;
  double fbest=0, fdifbest=0;
  int loopbest, itbest;
  arma::rowvec wdbest(K);
  
  bool args_ok = checkArgs(Q, Rndstart, verbose, maxiter, tol, J);
  bool k_ok = checkK(K, n);
  bool rot_ok = 1;
  if(rot != 0 && rot != 1){
    Rcpp::Rcout << "Error: rot must be an integer = 0 or = 1" << std::endl;
    rot_ok = 0;
  }
  
  Xs = preproc(Xs, prep);
  bool prep_ok = checkPrep(prep);
  bool stats_ok = checkStats(print);
  
  Xs2 = arma::pow(Xs,2);
  st = arma::sum(arma::sum(Xs2));
  
  S = arma::cov(Xs,1);
  S2x = arma::repmat(arma::sum(arma::pow(Xs,2),1),1,K); // assign function
  if(args_ok == 1 && k_ok ==1 && rot_ok ==1 && stats_ok && prep_ok == 1){
    for(loop=0;loop<Rndstart;loop++){
      // Initialization
      U = randPU(n, K);
      su = arma::sum(U, 0);
      Xmean = arma::diagmat(1.0/su)*U.t()*Xs;
      // Update A
      XX = Xs.t()*U*arma::diagmat(1.0/su)*U.t()*Xs;
      if(n>20000)
        arma::eig_sym(l, A, XX, "dc");
      else
        arma::eig_sym(l, A, XX);
      dL = arma::sort(l, "descend");
      idL = arma::sort_index(l, "descend");
      l = l(idL);
      l = l.subvec(0,Q-1);
      L = arma::diagmat(l);
      A = A.cols(idL);
      A = A.cols(0,Q-1);
      // Project Xs, Xmean on A
      Ymean = Xmean*A;
      Y = Xs*A;
      // f
      f0 = arma::trace(Ymean.t()*U.t()*U*Ymean);
      fdif=2*tol;
      it = 0;
      // Iteration phase
      while(fdif>tol || it >= maxiter){
        it++;
        S2x = arma::repmat(arma::sum(arma::pow(Y,2),1),1,K);
        U = assign_ssed(S2x, Ymean, Y, K, n);
        su = arma::sum(U,0);
        while(arma::sum(su==0)>0){ // If verified, solve the empty cluster issue
          Xmean = arma::diagmat(1.0/su)*U.t()*Xs;
          Xmean2 = arma::diagmat(1.0/su)*U.t()*Xs2;
          Ymean2 = arma::diagmat(1.0/su)*U.t()*arma::pow(Xs*A,2);
          wd = arma::trans(arma::sum((arma::diagmat(su)*(Ymean2 - arma::pow(Xmean*A,2))),1));
          wd.elem(arma::find(su==0 || su==1)).fill(0);
          U = split_maxwd(su, wd, U, Xs*A, K);
          su = arma::sum(U,0);
        }
        // Given U, compute Xmean (compute centroids)
        Xmean = arma::diagmat(1.0/su)*U.t()*Xs;
        // Given U and Xmean, update A
        XX = Xs.t()*U*arma::diagmat(1.0/su)*U.t()*Xs;
        
        if(n>20000)
          arma::eig_sym(l, A, XX, "dc");
        else
          arma::eig_sym(l, A, XX);
        dL = arma::sort(l, "descend");
        idL = arma::sort_index(l, "descend");
        l = l(idL);
        l = l.subvec(0,Q-1);
        L = arma::diagmat(l);
        A = A.cols(idL);
        A = A.cols(0,Q-1);
        // Project Xs and Xmean on A
        Ymean = Xmean*A;
        Y = Xs*A;
        // within sum of squares
        Ymean2 = arma::diagmat(1.0/su)*U.t()*arma::pow(Xs*A,2);
        Xmean2 = arma::diagmat(1.0/su)*U.t()*Xs2;
        wd = arma::trans(arma::sum((arma::diagmat(su)*(Ymean2 - arma::pow(Xmean*A,2))),1));
        // f
        f = arma::trace(Ymean.t()*U.t()*U*Ymean);
        fdif = f-f0;
        if(fdif>tol){
          f0=f;
          A0=A;
        }
        else{
          break;
        }
      }
      if(verbose==1)
        Rcpp::Rcout << "FKM: Loop = " << loop+1 << "; Explained variance (%) = " << (f/st)*100 << "; iter = " << it << "; fdif = " << fdif << std::endl;
      if(loop==0){
        Ubest=U;
        Abest=A;
        Ybest = Xs*Abest;
        Ymbest = Ymean;
        fbest = f;
        loopbest = 1;
        itbest = it+1;
        fdifbest = fdif;
        wdbest = wd;
      }
      if(f>fbest){
        Ubest=U;
        Abest = A;
        Ybest=Xs*Abest;
        Ymbest = Ymean;
        fbest=f;
        loopbest = loop+1;
        itbest = it+1;
        fdifbest = fdif;
        wdbest = wd;
      }
    }
    // sort components in descend order of variance
    // and rotate factors
    varY = arma::trans(var(Ybest,1));
    ic = arma::sort_index(varY, "descend");
    varY = varY(ic);
    Abest = Abest.cols(ic);
    Ybest = Ybest.cols(ic);
    //rotation of components axis
    if(rot != 0)
      Abest = varimax(Abest);
    // sort clusters of objects in descending order of cardinality
    iicc = arma::sort_index(arma::diagvec(Ubest.t()*Ubest), "descend");
    Ubest = Ubest.cols(iicc);
    wdbest = arma::trans(wdbest(iicc));
    double pseudoF = (f/(K-1))/((st-f)/(n-K));
    if(verbose==1)
      Rcpp::Rcout << "FKM (Final): Explained variance (%) = "<< (fbest/st)*100 <<"; loop = "<< loopbest << "; iter = " << itbest << "; fdif = " << fdifbest << "." << std::endl;
    if(print==1){
      Rcpp::Rcout << "\n>> Variance Explained by the FKM (% BSS / TSS): " << (fbest/st)*100 << std::endl;
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
      //for(int j = 0; j<J;j++){
      //  namesq(j) = "Variable" + Rcpp::toString(j+1);
      //}
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
      Rcpp::Rcout << "\n>> Summary of the latent factors:\n" << std::endl;
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
                            Rcpp::Named("centers") = Ymbest,
                            Rcpp::Named("withinss") = wdbest,
                            Rcpp::Named("betweenss") = fbest,
                            Rcpp::Named("totss") = st,
                            Rcpp::Named("size") = arma::sum(Ubest,0),
                            Rcpp::Named("pseudoF") = pseudoF,
                            Rcpp::Named("loop") = loopbest,
                            Rcpp::Named("it") = itbest);
  }
  else
    stop("Re-run with proper values for the arguments.");
}
