#' @name apseudoF
#' @title pseudoF (pF or Calinski-Harabsz) index for choosing k in partitioning models
#' @description
#' Calculates and plots the CH index for k = 2, ..., maxK. The function provides an interval wide (2tol*pF) so that the choice of K is less conservative. Instead of just choosing the maximum pF, if it exists, picks the value such that its upper bound is larger than max pF. 
#' 
#' @usage apseudoF(data, maxK, tol, model, Q)
#' 
#' @param data Units x variables numeric data matrix.
#' @param maxK Maximum number of clusters for the units to be tested.
#' @param tol Approximation value. It is half of the length of theinterval put for each pF. 0 <= tol < 1. Its default value is 0.05.
#' @param model Partitioning Models to run for each value of k. (1 = doublekm; 2 = redkm; 3 = factkm; 4 = dpcakm)
#' @param Q Number of principal components w.r.t. variables selected for the maxK -1 partitions to be tested.
#' 
#' @return \item{bestK}{best value of K (scalar).}
#' 
#' @author Ionel Prunila, Maurizio Vichi
#' 
#' @references 
#' Calinski T., Harabasz J. (1974) "A dendrite method for cluster analysis" <doi:10.1080/03610927408827101>
#' 
#' @examples
#' # Iris data 
#' # Loading the numeric variables of iris data
#' iris <- as.matrix(iris[,-5]) 
#' 
#' apF <- apseudoF(iris, maxK=10, tol = 0.05, model = 3, Q = 2)
#' 
#' @export

apseudoF = function(data, maxK=10, tol = 0.05, model, Q){
  if(any(missing(data), missing(model), missing(Q)))
    stop("The data set, model and Q must be given")
  if(any(is.null(data), is.null(model), is.null(Q)))
    stop("The data set, model and Q are empty")
  n=nrow(data)
  if(any(is.na(data), is.na(model), is.na(Q)))
    stop("The data set, model or Q must not contain NA values")
  if(any(!is.numeric(data), !is.numeric(model), !is.numeric(Q)))
    stop("The data set is not a numeric data.frame or matrix")
  if(maxK<2)
    stop("The maximum number of clusters to be tested must be >= 2")
  
  maxK <- as.integer(maxK)
  
  pseudoF <- vector(mode = "numeric", length = maxK-1)
  
  for(k in 2:maxK){
    if(model == 1)
      km <- doublekm(data, k, Q) 
    else if(model == 2 )
      km <- redkm(data, k, Q) 
    else if(model == 3)
      km <- factkm(data, k, Q)
    else if(model == 4)
      km <- dpcakm(data, k, Q)
    else
      stop("The value for model must be either 1, 2, 3 or 4.")
    pseudoF[k-1] <- km$pseudoF  
  }
  
  upseudoF = pseudoF*(1+tol)
  dpseudoF = pseudoF*(1-tol)
  im = which.max(pseudoF)
  maxpF = pseudoF[im]
  img <- im
  c <- 0
  for(i in im+1:length(pseudoF)){
    if(!is.na(upseudoF[i]) && !is.na(pseudoF[im]) && upseudoF[i]>=pseudoF[img]){
      img <- i
      c <- 1
      break
    }
  }
  if(c==1){
  im <- img}
  plot(2:maxK, pseudoF, type = "o", pch = NA, ylim = range(dpseudoF, upseudoF, pseudoF),
     xlab = "number of clusters", ylab = "", main = "Calinski-Harabasz Index (pseudoF)")

  segments(x0 = 2:maxK, y0 = pseudoF, x1 = 2:maxK, y1 = upseudoF, col = "blue")
  segments(x0 = 2:maxK, y0 = pseudoF, x1 = 2:maxK, y1 = dpseudoF, col = "red")
  message("The optimal number of clusters based on the pseudoF criterion is:", im+1, "\n")
  return("bestK" = im+1)
}