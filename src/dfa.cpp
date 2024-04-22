#define ARMA_WARN_LEVEL 0
#include "disfa.h"
#include "dispca.h"
#include "rc_aux.h"
#include "prep.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

//' @name disfa
//' @title Disjoint Factor Analysis
//' @description
//' Performs disjoint factor analysis, i.e., a Factor Analysis with a simple structure. In fact, each factor is defined by a disjoint subset of variables, resulting thus, in a simplified, easier to interpret loading matrix A and factors. Estimation is carried out via Maximum Likelihood.
//'  
//' 
//' @usage disfa(X, Q, Rndstart, verbose, maxiter, tol, constr, prep, print)
//' 
//' @param X Units x variables numeric data matrix.
//' @param Q Number of factors.
//' @param Rndstart Number of runs to be performed (Defaults is 20).
//' @param verbose Outputs basic summary statistics for each run (1 = enabled; 0 = disabled, default option).
//' @param maxiter Maximum number of iterations allowed (if convergence is not yet reached. Default is 100).
//' @param tol Tolerance threshold (maximum difference between the values of the objective function of two consecutive iterations such that convergence is assumed. Default is 1e-6).
//' @param constr is a vector of length J = nr. of variables, pre-specifying to which cluster some of the variables must be assigned. Each component of the vector can assume integer values from 1 o Q (See example for more details), or 0 if no constraint on the variable is imposed (i.e., it will be assigned based on the plain algorithm).
//' @param prep Pre-processing of the data. 1 performs the z-score transform (default choice); 2 performs the min-max transform; 0 leaves the data un-pre-processed.
//' @param print Prints summary statistics of the performed method (1 = enabled; 0 = disabled, default option).
//' 
//' @return returns a list of estimates and some descriptive quantities of the final results.
//' @return \item{V}{Variables x factors membership matrix (binary and row-stochastic). Each row is a dummy variable indicating to which cluster each variable has been assigned.}
//' @return \item{A}{Variables x components loading matrix.}
//' @return \item{Psi}{Specific variance of each observed variable, not accounted for by the common factors (matrix).}
//' @return \item{discrepancy}{Value of the objective function, to be minimized. Difference between the observed and estimated covariance matrices (scalar).}
//' @return \item{RMSEA}{Adjusted Root Mean Squared Error (scalar).}
//' @return \item{AIC}{Aikake Information Criterion (scalar).}
//' @return \item{BIC}{Bayesian Information Criterion (scalar).}
//' @return \item{GFI}{Goodness of Fit Index (scalar).}
//'
//' @author Ionel Prunila, Maurizio Vichi
//' 
//' @references 
//' Vichi M. (2017) "Disjoint factor analysis with cross-loadings" <doi:10.1007/s11634-016-0263-9>
//' 
//' 
//' @examples
//' # Iris data 
//' # Loading the numeric variables of iris data
//' iris <- as.matrix(iris[,-5]) 
//' 
//' # No constraint on variables
//' out <- disfa(iris, Q = 2)
//' 
//' # Constraint: the first two variables must contribute to the same factor.
//' outc <- disfa(iris, Q = 2, constr = c(1,1,0,0))
//' 
//' @export


//[[Rcpp::export]]

List disfa(Rcpp::NumericMatrix X ,int Q, int Rndstart = 10, int verbose = 0, int maxiter = 100, double tol = 1e-6, arma::vec constr = 00, int prep = 1, int print = 0){
  Rcpp::List Nomi = X.attr("dimnames");
  arma::mat Xs = as<arma::mat>(X);
  arma::mat S, VC = arma::eye(Q,Q), A, A0, AA, Psi, Sx, Sg, S0g, V, Vbest, Abest, Psibest, Ybest = {}, Y, L, Sxbest, invSxS, U;
  arma::vec a, aa, dAA, JCg, JCpm, dS, ibCg, ibCpm, varYbest, l, dL;
  arma::uvec indices, idL, ic;
  unsigned int J = Xs.n_cols;
  unsigned int n = Xs.n_rows;
  int loop, g, it;
  unsigned int j, posmin, flg, itbest=0, loopbest=0;
  double st, f0, fmin = arma::datum::inf, f=-arma::datum::inf, fdif=arma::datum::inf, ldS, GFI, AGFI, fbest, fdifbest = arma::datum::inf, df, llbest, pm, BIC, AIC, RMSEA, X2;
  
  bool args_ok = checkArgs(Q, Rndstart, verbose, maxiter, tol, J);
  bool constr_ok = checkConstr(constr, J, Q);
  Xs = preproc(Xs, prep);
  bool prep_ok = checkPrep(prep);
  bool stats_ok = checkStats(print);
  if(constr.size()==1)
    constr = arma::zeros(J,1);
  if(args_ok ==1 && constr_ok == 1 && prep_ok ==1 && stats_ok ==1){
    //var-covar matrix
    if(J!=n)
      S = arma::cov(Xs, 1);
    else
      S = Xs;
    dS = arma::diagvec(S);
    ldS = log(arma::det(S));
    if(J!=n){
      if(arma::rank(Xs)<J){
        Rcpp::Rcout << "Attention! Singular Variance-Covariance Matrix." << std::endl;
        ldS = 1;
      }
      st = (1.0/n)*arma::sum(arma::sum(arma::pow(Xs,2)));//
    }
    else
      st = arma::trace(Xs);
    arma::vec JJ = arma::regspace(0,1,J-1);
    arma::vec zJ = arma::zeros(J,1);
    //start
    for(loop=0; loop < Rndstart; loop++){
      it=0;
      flg = 1;
      while(flg>0){
        V = randPU(J,Q);
        for(j=0; j < J; j++){
          if(constr(j)>0)
            V.row(j) = VC.row(constr(j)-1);
        }
        flg = arma::sum(arma::find(arma::sum(V,0)==0));
      }
      A = arma::zeros(J,Q);
      
      for(g=0; g<Q; g++){
        ibCg = V.col(g);
        JCg = JJ.elem(arma::find(ibCg==1));
        if(sum(ibCg)>1){
          Sg = S.submat(arma::conv_to<arma::uvec>::from(JCg), arma::conv_to<arma::uvec>::from(JCg));
          a = AF1(Sg);
          A.submat(arma::conv_to<arma::uvec>::from(JCg),arma::uvec{static_cast<unsigned int>(g)}) = a;
        }
        else{
          A.submat(arma::conv_to<arma::uvec>::from(JCg), arma::uvec{static_cast<unsigned int>(g)}) = arma::ones(JCg.n_elem);
        }
      }
      
      A0 = A;
      // initial Psi
      AA = A*A.t();
      dAA = arma::diagvec(AA);
      Psi = arma::diagmat(dS-dAA);
      Sx = AA + Psi;
      f0 = log(arma::det(Sx)) - ldS + arma::trace(solve(Sx,S, arma::solve_opts::allow_ugly)) - J;
      fmin = arma::datum::inf;
      // iteration phase
      for(it=0; it<maxiter; it++){
        // update V and A
        for(j=0; j<J; j++){
          posmin = JJ(arma::index_max(V.row(j)==1));
          if(constr(j)==0){
            for(g = 0; g < Q; g++){
              V.row(j) = VC.row(g);
              if(arma::sum(V.col(posmin))>0){
                ibCg = V.col(g); // new class of V
                ibCpm = V.col(posmin); // old class of V
                indices = arma::find(ibCg==1);
                JCg = JJ.elem(indices);
                indices = arma::find(ibCpm==1);
                JCpm = JJ.elem(indices);
                
                
                if(sum(ibCg)>1){
                  Sg = S.submat(arma::conv_to<arma::uvec>::from(JCg), arma::conv_to<arma::uvec>::from(JCg));
                  a = AF1(Sg);
                  if(arma::sum(a)<0){
                    a = -a;
                  }
                  A.col(g) = zJ;
                  A.submat(arma::conv_to<arma::uvec>::from(JCg), arma::uvec{static_cast<unsigned int>(g)}) = a;
                }
                else{
                  A.col(g) = zJ;
                  A.submat(arma::conv_to<arma::uvec>::from(JCg), arma::uvec{static_cast<unsigned int>(g)}) = arma::ones(JCg.n_elem);
                }
                if(sum(ibCpm)>1){
                  S0g = S.submat(arma::conv_to<arma::uvec>::from(JCpm), arma::conv_to<arma::uvec>::from(JCpm));
                  aa = AF1(S0g);
                  if(arma::sum(aa)<0){
                    aa = -aa;
                  }
                  A.col(posmin) = zJ;
                  A.submat(arma::conv_to<arma::uvec>::from(JCpm), arma::uvec{static_cast<unsigned int>(posmin)}) = aa;
                }
                else{
                  A.col(posmin) = zJ;
                  A.submat(arma::conv_to<arma::uvec>::from(JCpm), arma::uvec{static_cast<unsigned int>(posmin)}) = arma::ones(JCpm.n_elem);
                }
                AA = A*A.t();
                dAA = arma::diagvec(AA);
                Psi = arma::diagmat(dS-dAA);
                Sx = AA + Psi;
                f = log(arma::det(Sx)) - ldS + arma::trace(solve(Sx,S, arma::solve_opts::allow_ugly)) - J;//
                if(f < fmin){
                  fmin = f;
                  posmin = g;
                  A0 = A;
                }
                else{
                  A = A0;
                }
              }
            }
          }
          V.row(j) = VC.row(posmin);
        }
        
        AA = A*A.t();
        dAA = arma::diagvec(AA);
        Psi = arma::diagmat(dS-dAA);
        Sx = AA + Psi;
        f = log(arma::det(Sx)) - ldS + arma::trace(solve(Sx,S, arma::solve_opts::allow_ugly)) - J;
        fdif = f0-f;
        if(fdif > tol){
          f0=f;
          fmin=f0;
          A0=A;
        }
        else{
          break;
        }
      }
      invSxS = solve(Sx,S, arma::solve_opts::allow_ugly); 
      GFI = 1 - trace((invSxS-arma::eye(J,J))*(invSxS-arma::eye(J,J)))*(1.0/(arma::trace((invSxS)*(invSxS))));
      df = (pow((J-Q),2) + (J+Q))*0.5; 
      AGFI = 1 - (J*(J+1)/(2*df))*(1-GFI);
      //
      if(verbose==1)
        Rcpp::Rcout << "DFA: Loop = " << loop+1 << "; Discrepancy = " << f << "; Goodness of Fit = " << GFI << "; iter = "<< it+1 << "; fdif = " << fdif << std::endl;
      if(loop==0){
        Vbest = V;
        Abest = A;
        Psibest = arma::diagvec(Psi);
        if(J!=n){
          Ybest = zscore(Xs*Abest)*arma::pow(Abest.t()*Abest, 0.5);
        }
        else{
        }
        fbest = f;
        loopbest = 1;
        itbest = it+1;
        fdifbest = fdif;
      }
      if(f < fbest){
        Vbest = V;
        fbest = f;
        Abest = A;
        Psibest = arma::diagvec(Psi);
        if(J!=n){
          Ybest = zscore(Xs*Abest)*arma::pow(Abest.t()*Abest, 0.5);
        }
        else{
        }
        loopbest = loop+1;
        itbest = it+1;
        fdifbest = fdif;
      }
      else if(f-fbest<0.0001 && loop > 0){
        //break;
      }
    }
    if(sum(sum(Abest,0)==0)>0)
      Rcpp::Rcout << "There is an empty cluster of variables. Consider setting to zero more variables in the constr argument." << std::endl;
    
    // sort the final solution in descending order of variance
    if(J!=n){
      varYbest = arma::trans(arma::var(Ybest,1));
      ic = arma::sort_index(varYbest, "descend");
      Abest = Abest.cols(ic);
      Vbest = Vbest.cols(ic);
      Ybest = Ybest.cols(ic);
    }
    // variance of the second component
    arma::vec e2k = arma::zeros(Q,1);
    arma::vec cro = arma::zeros(Q,1);
    for(g = 0; g<Q; g++){
      ibCg = Vbest.col(g);
      if(arma::sum(ibCg)>1){
        indices = arma::find(ibCg==1);
        JCg = JJ.elem(indices);
        cro(g) = CronbachAlpha(Xs.cols(arma::conv_to<arma::uvec>::from(JCg)));
        Sg = S.submat(arma::conv_to<arma::uvec>::from(JCg), arma::conv_to<arma::uvec>::from(JCg));
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
    //compute
    Sxbest = Abest*Abest.t() + arma::diagmat(Psibest);
    invSxS = solve(Sx,S, arma::solve_opts::allow_ugly);
    GFI = 1 - arma::trace((invSxS-arma::eye(J,J))*(invSxS-arma::eye(J,J)))*(1.0/(arma::trace((invSxS)*(invSxS))));
    llbest = -(0.5)*n*log(2*arma::datum::pi)-(0.5)*n*(log(arma::det(Sxbest)) + arma::trace(invSxS));
    pm = 2*J-Q;
    df = J*(J-3)/2+Q;
    AGFI = 1 - (J*(J+1)/(2*df))*(1-GFI);
    X2 = (n-(2*J+4*Q+11)/6)*log(arma::det(Sxbest)*(1.0/arma::det(S)));
    //pValue
    BIC = -2*llbest + log(n) * pm;
    AIC = -2*llbest + 2*pm;
    arma::vec t(2);
    t(0) = 0;
    t(1) = fbest/(df)-1/(n-1);
    RMSEA = std::sqrt(max(t));
    double se;
    arma::vec p;
    if(verbose==1)
      Rcpp::Rcout << "DFA (Final): Loop " << loopbest << "; Discrepancy = " << fbest << "; Goodness of Fit = " << GFI <<"; Adjusted Goodness of Fit = " << AGFI <<"; iter = "<< itbest << "; fdif = " << fdifbest << "\n"<<std::endl;
      //Rcpp::Rcout << << std::endl;
    if(print == 1){
      Rcpp::Rcout << ">> Discrepancy of DFA: " << fbest << std::endl;
      Rcpp::Rcout << "\n>> Summary statistics:\n" << std::endl;
      Rcpp::CharacterMatrix Summary(1,6);
      Rcpp::NumericVector valori = {pm, X2, df, BIC, AIC, RMSEA};
      Summary(0,_) = valori;
      Rcpp::CharacterVector names = {"Unknown Parameters", "Chi-square", "Degrees of Freedom", "BIC", "AIC", "RMSEA"};
      Rcpp::colnames(Summary) = names;
      names = " ";
      Rcpp::rownames(Summary) = names;
      Function print("print");
      print(Summary, Named("quote") = false);
      Rcpp::CharacterVector vars = Nomi[1]; 
      Rcpp::Rcout << "\n>> Loading Matrix (Manifest Variables x Latent Variables) \n" << std::endl;
      Rcpp::NumericMatrix matrice = wrap(Abest);
      Rcpp::rownames(matrice) = vars;
      Rcpp::CharacterVector facts(Q);
      for(int j=0; j<Q; j++){
        facts(j) = "Factor " + Rcpp::toString(j+1);
      }
      Rcpp::colnames(matrice) = facts;
      Rf_eval(Rf_lang2(Rf_install("print"), matrice), R_GlobalEnv);
      double varcum = 0;
      arma::vec com;
      Rcpp::Rcout << "\n>> Summary of the latent factors:\n" << std::endl;
      Rcpp::CharacterVector fields(6);
      fields = {"Explained Variance", "Expl. Var. (%)", "Cum. Var", "Cum. Var (%)", "Var. 2nd component", "Cronbach's Alpha"};
      Rcpp::NumericMatrix fattori(Q,6);
      for(int q = 0; q<Q; q++){
        varcum = varcum + arma::as_scalar(sum(pow(Abest.col(q),2)));
        valori = {arma::as_scalar(sum(pow(Abest.col(q),2))), (arma::as_scalar(sum(pow(Abest.col(q),2)))/st)*100, varcum,( (varcum/st)*100), e2k(q), cro(q)};
        fattori(q,_) = valori;
      }
      Rcpp::colnames(fattori) = fields;
      Rcpp::rownames(fattori) = facts;
      Rf_eval(Rf_lang2(Rf_install("print"), fattori), R_GlobalEnv);
      com = sum(pow(Abest,2),1);
      Rcpp::Rcout << "\n>> Detailed Manifest-variable - Latent-factor relationships\n" << std::endl;
      fields = {"Associated Factor", "Corr. Coeff.", "Std. Error", "Pr(p>|Z|)", "Var. Error", "Communality"};
      Rcpp::NumericMatrix variabili(J,6);
      for(j = 0; j < J; j++){
        p = arma::zeros(J,1);
        for(double q = 0; q<Q; q++){
          if(Abest(j,q)!=0){
            p(j) = r2pv(Abest(j,q),n);
            se = std::sqrt(std::abs(Psibest(j))) * (1/std::sqrt(n));
            valori = {q+1, Abest(j,q), se, p(j), Psibest(j), com(j)};
          }
          variabili(j,_) = valori;
        }
      }
      Rcpp::colnames(variabili) = fields;
      Rcpp::rownames(variabili) = vars;
      Rf_eval(Rf_lang2(Rf_install("print"), variabili), R_GlobalEnv);
    }
  return Rcpp::List::create(Rcpp::Named("GFI") = GFI,
                            Rcpp::Named("AGFI") = AGFI,
                            Rcpp::Named("A") = Abest,
                            Rcpp::Named("V") = Vbest,
                            Rcpp::Named("Psi") = Psibest,
                            Rcpp::Named("discrepancy") = fbest,
                            Rcpp::Named("RMSEA") = RMSEA,
                            Rcpp::Named("AIC") = AIC,
                            Rcpp::Named("BIC") = BIC,
                            Rcpp::Named("X2") = X2);
  }
  else
    stop("Re-run with proper values for the arguments.");
}
