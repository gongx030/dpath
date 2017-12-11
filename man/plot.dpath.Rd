\name{plot.dpath}
\alias{plot.dpath}
\title{
Visualizing a \code{dpath} object
}
\description{
Visualizing a \code{dpath} object.
}
\usage{
\method{plot}{dpath}(x, ...)
}
\arguments{
  \item{x}{a \code{dpath} object}
	\item{\dots}{Arguments to be passed to methods:
		\code{type} what type of plot should be drawn.  Possible types are
		\itemize{
			\item \code{markers}: metagene coefficients, metagene basis
			and observed expression levels for selected marker genes.
			\item \code{cell.cluster}: results of cell clustering by partioning the
			metacell landscape, as well as the metagene coefficients for each cluster.
			\item \code{metagene.entropy}: metagene entropy on the SOM.
			\item \code{metacell.landscape}: the metacell landscape, with
			user-specified properties for highlighting different metacells
			\item \code{prioritization}: prioritization score of metacells with respect 
			to specified cellular state and top ranking genes which expression
			pattern is most correlated with the prioritization score
			\item \code{gene.expression}: gene expression pattern on the SOM.
		}
		\code{type} specific arguments:
		\itemize{
		  \item \code{genes}: gene name(s); required when \code{type} is \code{markers} or
			\code{gene.expression}
			\item \code{cell.group}: cell groups; required when
			\code{type} is \code{cell.distribution} or \code{metagene.entropy}
			\item \code{score.gene}: gene ranking; required when \code{type} is
			\code{prioritization}
			\item \code{score.metacell} cell ranking; required when \code{type} is
			\code{prioritization}
		}
	}
}
\author{
Wuming Gong <gongx030@umn.edu>
}
\examples{
	help('dpath')
}
\keyword{dpath}
