\name{fitsom}
\alias{fitsom}
\title{
Fitting a self-organizing map (SOM) using bootstrapped cells and clustering
metacells by partitinoning the metacell landscape
}
\description{
This function first build a self organizing map (SOM) as a metacell landscape using the metagene
coefficients of bootstrapped cells, and partion the metacell landscape into
smaller regions.
}
\usage{
fitsom(dp, k = NA, xdim = 15, ydim = 15, n.min = 10)
}
\arguments{
  \item{dp}{A \code{dpath} object}
	\item{k}{The desired number of clusters}
	\item{xdim,ydim}{Dimensions of the metacell landscape}
	\item{n.min}{The minimum number of cells in a cluster}
}
\value{
	A \code{dpath} object
}
\author{
Wuming Gong <gongx030@umn.edu>
}
\examples{
	help('dpath')
}
\keyword{dpath}
