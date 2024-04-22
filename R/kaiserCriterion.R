#' @name kaiserCrit
#' @title Selecting the number of principal components to be extracted from a dataset
#' @description
#' Selects the optimal number of principal components to be extracted from a dataset based on Kaiser's criterion
#' 
#' @usage kaiserCrit(data)
#' 
#' @param data Units x variables data matrix.
#' 
#' @return \item{bestQ}{Number of components to be extracted (scalar).}
#' 
#' @author Ionel Prunila, Maurizio Vichi
#' 
#' @references 
#' Kaiser H. F. (1960) "The Application of Electronic Computers to Factor Analysis" <doi:10.1177/001316446002000>
#' 
#' @examples
#' # Iris data 
#' # Loading the numeric variables of iris data
#' iris <- scale(as.matrix(iris[,-5])) 
#' 
#' # Apply the Kaiser rule
#' h <- kaiserCrit(iris)
#' 
#' @export

kaiserCrit = function(data){
  data <- as.matrix(data)
  if(missing(data))
    stop("The dataset must be given")
  if(is.null(data))
    stop("The dataset is empty")
  n=nrow(data)
  if(!is.numeric(data))
    stop("The data set is not a numeric data.frame or matrix")
  
  data <- scale(data)
  eigs <- eigen(cov(data))
  Q = sum(eigs$values>0.9)
  message("The number of components suggested by the Kaiser criterion is: ", Q, "\n")
  return("bestQ" = Q)
}
