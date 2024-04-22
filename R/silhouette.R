#' @name silhouette
#' @title Silhouette
#' @description
#' Computes and plots the silhouette of a partition
#' 
#' @usage silhouette(data, drclust_out)
#' 
#' @param data Units x variables data matrix.
#' @param drclust_out Out of either doublekm, redkm, factkm or dpcakm.
#' 
#' @return \item{cl.silhouette}{Silhouette index for the given partition, for each object (matrix).}
#' @return \item{fe.silhouette}{Factoextra silhouette graphical object}
#' 
#' @import cluster
#' @import factoextra
#' @author Ionel Prunila, Maurizio Vichi
#' 
#' @references 
#' Rousseeuw P. J. (1987) "Silhouettes: a Graphical Aid to the Interpretation and Validation of Cluster Analysis" <doi:10.1016/0377-0427(87)90125-7>
#' 
#' Maechler M. et al. (2023) "cluster: Cluster Analysis Basics and Extensions" <https://CRAN.R-project.org/package=cluster>
#' 
#' Kassambara A. (2022) "factoextra: Extract and Visualize the Results of Multivariate Data Analyses" <https://cran.r-project.org/web/packages/factoextra/index.html>
#' 
#' @examples
#' # Iris data 
#' # Loading the numeric variables of iris data
#' iris <- as.matrix(iris[,-5]) 
#' 
#' #standardizing the data
#' iris <- scale(iris)
#' 
#' #applying a clustering algorithm
#' drclust_out <- dpcakm(iris, 20, 3)
#' 
#' #silhouette based on the data and the output of the clustering algorithm
#' d <- silhouette(iris, drclust_out)
#' 
#' @export

require(cluster)
require(factoextra)
silhouette = function(data, drclust_out){
  if(any(missing(data),missing(drclust_out)))
    stop("The dataset and the model must be given")
  if(any(is.null(data), is.null(drclust_out)))
    stop("The data or the model are empty")
  n=nrow(data)
  if(any(is.na(data), is.na(drclust_out)))
    stop("The data and the model must not contain NA values")
  if(any(!is.numeric(data)))
    stop("The data set is not a numeric data.frame or matrix")
  
  sil <- cluster::silhouette(cluster(drclust_out$U), dist(data)) 
  fsil <- factoextra::fviz_silhouette(sil)
  return(list("cl.silhouette" = sil, "fe.silhouette" = fsil))
}