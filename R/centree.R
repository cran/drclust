#' @name centree
#' @title Ward-dendrogeam of centroids of partitioning models
#' @description
#' Plots the Ward-dendrogram of the centroids of a partitioning model. The plot is useful as a diagnosis tool for the choice o the number of clusters.
#' 
#' @usage centree(drclust_out)
#' 
#' @param drclust_out Output of either doublekm, redkm, factkm or dpcakm.
#' 
#' @return \item{centroids-dkm}{Centroids x centroids distance matrix.}
#' 
#' @author Ionel Prunila, Maurizio Vichi
#' 
#' @references 
#' Ward J. H. (1963) "Hierarchical Grouping to Optimize an Objective Function" <doi:10.1080/01621459.1963.10500845>
#' 
#' @examples
#' # Iris data 
#' # Loading the numeric variables of iris data
#' iris <- as.matrix(iris[,-5]) 
#' 
#' dc_out <- dpcakm(iris, 20, 3)
#' d <- centree(dc_out)
#' 
#' @export

centree = function(drclust_out){
  tryCatch(
    expr = {
  centers = drclust_out$centers
  if(nrow(centers)<3)
    stop("The dendrogram can be built from at least 3 centroids")
  if(missing(centers))
    stop("The centroids must be given")
  if(is.null(centers))
    stop("The centroid matrix is empty")
  if(any(is.na(centers)))
    stop("The centroid matrix must not contain NA values")
  if(!is.numeric(centers))
    stop("The centroid matrix is not numeric")
  d <- dist(centers)
  hc <- hclust(d, method = "ward.D")
  plot(hc, xlab = "centroids", ylab = "Increase in the loss-function", main = "Centroids Dendrogram")
  return("centroids-dm" = d)
    }
  ,
  warning = function(w){
    message("Please provide the proper argument to the function")
  }
  )
}
