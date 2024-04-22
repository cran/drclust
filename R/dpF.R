#' @name dpseudoF
#' @title double pseudoF (Calinski-Harabsz) index
#' @description
#' A pseudoF version for double partitioning, for the choice of the number of clusters of the units and variables (rows and columns of the data matrix). It is a diagnostic tool for inspecting simultaneously the optimal number of unit-clusters and variable-clusters.
#' 
#' @usage dpseudoF(data, maxK, maxQ)
#' 
#' @param data Units x variables numeric data matrix.
#' @param maxK Maximum number of clusters for the units to be tested.
#' @param maxQ Maximum number of clusters for the variables to be tested.
#' 
#' @return \item{dpseudoF}{matrix containing the pF value for each pair of K and Q within the specified range}
#' 
#' @author Ionel Prunila, Maurizio Vichi
#' 
#' @references 
#' R. Rocci, M. Vichi (2008)" Two-mode multi-partitioning" <doi:10.1016/j.csda.2007.06.025>
#' 
#' 
#' T. Calinski & J. Harabasz (1974). A dendrite method for cluster analysis. Communications in Statistics, 3:1, 1-27
#' 
#' @examples
#' # Iris data 
#' # Loading the numeric variables of iris data
#' iris <- as.matrix(iris[,-5]) 
#' 
#' dpeudoF <- dpseudoF(iris, maxK=10, maxQ = 3)
#' 
#' @export

dpseudoF = function(data, maxK=10, maxQ=10){
  if(any(missing(data), missing(maxQ), missing(maxK)))
    stop("The data set, maxK and maxQ must be given")
  if(any(is.null(data), is.null(maxK), is.null(maxQ)))
    stop("The data set, maxK and maxQ are empty")
  n=nrow(data)
  if(any(is.na(data), is.na(maxQ), is.na(maxK)))
    stop("The data set, maxK or maxQ must not contain NA values")
  if(any(!is.numeric(data), !is.numeric(maxQ), !is.numeric(maxK)))
    stop("The data set, maxQ o maxK is not a numeric data.frame or matrix")
  if(maxK<2 && maxK< nrow(data))
    stop("The maximum number of unit-clusters to be tested must be >= 2 and <= nrow(data)")
  if(maxQ<2 && maxQ< ncol(data))
    stop("The maximum number of variable-clusters to be tested must be >= 2 and <= ncol(data)")
  maxK <- as.integer(maxK)
  pseudoF <- matrix(nrow = maxK-1, ncol = maxQ-1)
  J = ncol(data)
  for(k in 2:maxK){
    for(q in 2:maxQ){
      dkm = doublekm(data, k, q)
      pseudoF[k-1,q-1] = (norm((dkm$U%*%solve(t(dkm$U)%*%dkm$U)%*%t(dkm$U)%*%scale(data)%*%dkm$V%*%solve(t(dkm$V)%*%dkm$V)%*%t(dkm$V) - (1/(n*J))*matrix(1, ncol = n, nrow = n)%*%scale(data)%*% matrix(1, ncol = J, nrow = J)), type = "F")/((k*q)-1))/(norm((dkm$U%*%solve(t(dkm$U)%*%dkm$U)%*%t(dkm$U)%*%scale(data)%*%dkm$V%*%solve(t(dkm$V)%*%dkm$V)%*%t(dkm$V)),type = "F")/(n*J - k*q))
    }
  }
  namesr = vector(mode = "character", length = maxK-1)
  namesc = vector(mode = "character", length = maxQ-1)
  for(k in 2:maxK){
    namesr[k-1] = paste("K =", k)
    #rownames(pseudoF)[k-1] = paste("K = ", k)
  }
  for(q in 2:maxQ){
    namesc[q-1] = paste("Q =", q)
    #colnames(pseudoF)[q-1] = paste("Q = ", q)
  }
  colnames(pseudoF) = namesc
  rownames(pseudoF) = namesr
  return(pseudoF)
}