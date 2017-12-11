\name{differentiation.path.p2c}
\alias{differentiation.path.p2c}
\title{
Find the shortest differentiation paths from predicted progenitor states to committed states
}
\description{
This function first predicts the most likely progenitor and committed states with respect to the 
user specified metagenes and find the shortest paths between them on the metacell landscape
}
\usage{
differentiation.path.p2c(dp, metagene = NULL)
}
\arguments{
  \item{dp}{A \code{dpath} object}
	\item{metagene}{The index of specified metagenes}
}
\value{
	A list: 
	\itemize{
		\item \code{mcp} The index of the most likely progenitor metacell
		\item \code{paths} A list of vectors indicating the paths from mcp to most likely committed metacells
		\item \code{A} An igraph object for resulting differentiation paths
	}
}
\author{
Wuming Gong <gongx030@umn.edu>
}
\examples{
	help('dpath')
}
\keyword{dpath}
