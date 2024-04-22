# data heatmap on a reduced subspace (ordered by clusters, which are also ordered)

#' @name heatm
#' @title Heatmap of a partition in a reduced subspace
#' @description
#' Plots the heatmap of a partition on a reduced subspace obtained via either: doublekm, redkm, factkm or dpcakm.
#' 
#' @usage heatm(data, drclust_out)
#' 
#' @param data Units x variables data matrix.
#' @param drclust_out Out of either doublekm, redkm, factkm or dpcakm.
#' 
#' @return No return value, called for side effects
#'
#' @author Ionel Prunila, Maurizio Vichi
#' 
#' @references 
#' Kolde R. (2019) "pheatmap: Pretty Heatmaps" <https://cran.r-project.org/web/packages/pheatmap/index.html> 
#' 
#' @import pheatmap
#' @examples
#' # Iris data 
#' # Loading the numeric variables of iris data
#' iris <- as.matrix(iris[,-5]) 
#' 
#' # standardizing the data
#' iris <- scale(iris)
#' 
#' # applying a clustering algorithm
#' drclust_out <- dpcakm(iris, 20, 3)
#' 
#' # obtain a heatmap based on the output of the clustering algorithm and the data
#' h <- heatm(iris, drclust_out)
#' 
#' @export

require(pheatmap)
heatm = function(data, drclust_out){
  
  if(any(missing(data),missing(drclust_out)))
    stop("The dataset and the model must be given")
  if(any(is.null(data), is.null(drclust_out)))
    stop("The data or the model are empty")
  n=nrow(data)
  if(any(is.na(data), is.na(drclust_out)))
    stop("The data and the model must not contain NA values")
  if(any(!is.numeric(data)))
    stop("The data set is not a numeric data.frame or matrix")
  
  ord_data <- matrix(0, ncol = ncol(data))
  ord_data <- ord_data[-1,]
  
  K <- nrow(drclust_out[["centers"]])
  
  oDU <- matrix(0,ncol = K)
  oDU <- oDU[-1,]
  n = nrow(data%*%drclust_out[[2]])
  K = nrow(drclust_out[["centers"]])
  DU = matrix(0,nrow = n, ncol = K)
  for(i in 1:n){
    for(j in 1:K){
      if(drclust_out[["U"]][i,j]==1)
        DU[i,j] = sum((data[i,]%*%drclust_out[[2]] - drclust_out[["centers"]][j,])^2)
      #else
      #  DU[i,j]==0
    }
  }
  #DU <- distm(data%*%drclust_out[[2]], drclust_out[["centers"]])*drclust_out[["U"]]
  #D <- distm(data%*%drclust_out[[2]], drclust_out[["centers"]])
  ds <- sort(rowSums(drclust_out[["centers"]]), index.return = T)$ix

  for(i in ds){
    ind <- which(DU[,i]!=0)
    sind <- sort(DU[ind,i], index.return = T)$ix

    oDU <- rbind(oDU, DU[ind[sind],])
    ord_data <- rbind(ord_data, data[ind[sind],])
  }
  pheatmap::pheatmap(ord_data%*%drclust_out[[2]], cexCol = 0.7, main = "Heatmap", cluster_cols = FALSE)
}
