\name{prioritize}
\alias{prioritize}
\title{
Prioritizing metacells with respect to spceified cellular states and ranking genes by enrichment scores
}
\description{
Prioritize the metacells with respect to specified cellular states (progenitor
or committed) and rank genes based on the correlation between their expression
pattern and cellular states of the metacells.  
}
\usage{
prioritize(dp, direction = NULL, metagene = NULL, beta = 0.9)
}
\arguments{
  \item{dp}{A \code{dpath} object}
  \item{direction}{Cellular state:
		\itemize{
			\item \code{progenitor} progenitor or multi-potent state
			\item \code{committed} committed or differentiated state
		}
	}
	\item{metagene}{The desired weight for each metagene} 
	\item{beta}{Probability of staying in the metagene graph}
}
\value{
	A list of two vectors:
	\itemize{
		\item \code{gene} gene ranking score
		\item \code{metacell} metacell prioritization score
	}
}
\author{
Wuming Gong <gongx030@umn.edu>
}
\examples{
	help('dpath')
}
\keyword{dpath}
