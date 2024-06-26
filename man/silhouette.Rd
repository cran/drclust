% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/silhouette.R
\name{silhouette}
\alias{silhouette}
\title{Silhouette}
\usage{
silhouette(data, drclust_out)
}
\arguments{
\item{data}{Units x variables data matrix.}

\item{drclust_out}{Out of either doublekm, redkm, factkm or dpcakm.}
}
\value{
\item{cl.silhouette}{Silhouette index for the given partition, for each object (matrix).}

\item{fe.silhouette}{Factoextra silhouette graphical object}
}
\description{
Computes and plots the silhouette of a partition
}
\examples{
# Iris data 
# Loading the numeric variables of iris data
iris <- as.matrix(iris[,-5]) 

#standardizing the data
iris <- scale(iris)

#applying a clustering algorithm
drclust_out <- dpcakm(iris, 20, 3)

#silhouette based on the data and the output of the clustering algorithm
d <- silhouette(iris, drclust_out)

}
\references{
Rousseeuw P. J. (1987) "Silhouettes: a Graphical Aid to the Interpretation and Validation of Cluster Analysis" <doi:10.1016/0377-0427(87)90125-7>

Maechler M. et al. (2023) "cluster: Cluster Analysis Basics and Extensions" <https://CRAN.R-project.org/package=cluster>

Kassambara A. (2022) "factoextra: Extract and Visualize the Results of Multivariate Data Analyses" <https://cran.r-project.org/web/packages/factoextra/index.html>
}
\author{
Ionel Prunila, Maurizio Vichi
}
