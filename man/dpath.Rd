\name{dpath}
\alias{dpath}
\title{
Decomposing expression matrix into metagenes and metacells
}
\description{
  This function takes the expression profile matrix of single cell RNA-seq as
  the input decomposes it into metagenes by Poisson mixture non-negative matrix
  factorization with the awareness of the dropout events, and cluster cells into
  metacells by metagene expression profiles. 
}
\usage{
dpath(X, K = 5, subset.gene = NULL, subset.cell = NULL, lambda0 = 0.1, 
max.iter = 500, repeat.mf = 50, mc.cores = 8)
}
\arguments{
  \item{X}{
  An expression profile matrix where each row represents a gene and each column
  represents a cell. The expression data need to be represented as 
	transcript per million (TPM). 
}
  \item{K}{Number of desired metagenes}
  \item{subset.gene}{a subset of genes for initializing factorization}
  \item{subset.cell}{a subset of cell for initializing factorization}
  \item{lambda0}{Mean expression of the dropout events}
  \item{max.iter}{Maximum rounds of iterations}
  \item{repeat.mf}{Repetitive runs of weighted Poisson NMF}
  \item{mc.cores}{Number of utilized CPU cores}
}
\value{
	A \code{dpath} object
}
\references{
	"Dpath software reveals hierarchical hemato-endothelial lineages of Etv2-
	progenitors based on single cell transcriptome analysis", W Gong, TL Rasmussen, 
	N Kayano-Nakagawa, BN Singh, W Pan and DJ Garry, 2016
}
\author{
Wuming Gong <gongx030@umn.edu>
}
\examples{
library(dpath)

# loading an expression profile matrix of 548 genes and 100 cells
data(etv2)
cell.group <- factor(colnames(etv2), c('E7.25', 'E7.75', 'E8.25'))

# run wp-NMF with 4 metagenes, using cells from E7.75 and E8.25 to initialize
# the factorization, repeat the factorization for 50 times and utilize 8 CPU
# cores
# dp <- dpath(etv2, K = 4, subset.cell = cell.group %in% c('E7.75', 'E8.25'), 
# repeat.mf = 50, mc.cores = 8)

# reorder the metagenes from different repetitive runs
#dp <- reorder(dp)

# load the finished dp file
data(dp)

# visualizing metagene coefficients, metagene basis and observed expression
# levels for selected marker genes
gene.list <- c('Gata1', 'Ikzf1', 'Itga2b', 'Hba-a1', 'Runx1', 'Gata4', 
'Smarcd3', 'Tbx20', 'Alcam', 'Cgnl1', 'Dok4', 'Plvap', 'Emcn', 'Pecam1', 
'Cd34', 'Cdh5', 'Icam1', 'T', 'Kdr', 'Pdgfra', 'Gli2', 'Pou5f1', 'Nanog')	
dev.new(width = 22, height = 15)
plot(dp, type = 'markers', genes = gene.list, reorder.genes = FALSE)

# fitting a self-organizing map (SOM) using bootstrapped cells and clustering
# metacells by partitinoning the metacell landscape
set.seed(6580)
dp <- fitsom(dp, xdim = 15, ydim = 15, n.min = 15)

# plot the cell clustering results and the mean metagene coefficients of each
# cluster
dev.new(width = 9, height = 6)
plot(dp, type = 'cell.cluster', cell.group = cell.group)

# visualizing metagene entropy on the SOM
dev.new(width = 15, height = 15)
plot(dp, type = 'metagene.entropy', cell.group = cell.group)

# visualizing expression pattern of selected genes on the SOM
dev.new(width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(3, 3, 5, 1))
plot(dp, type = 'gene.expression', genes = c('Runx1', 'Gata1', 'Etv2', 
'Plvap', 'T', 'Tbx20'))

# prioritizing metacells with respect to the committed state of the 1st metagene
# (endothelial metagene)
pg.endothelial <- prioritize(dp, metagene = c(1, 0, 0, 0), direction = 'committed')
library(gplots)
col.metacell <- colorpanel(100, low = 'black', mid = 'white', high = 'purple')

# visualizing the prioritization score of committed endothelial lineages (the
# 1st metagene)
dev.new(width = 7, height = 7)
plot(dp, type = 'metacell.landscape', property = pg.endothelial$metacell, 
col.metacell = col.metacell)

# show the top 100 most enriched genes in the committed endothelial lineage (the 1st
# metagene)
dev.new(width = 22, height = 15)
plot(dp, type = 'prioritization', score.gene = pg.endothelial[['gene']], 
score.metacell = pg.endothelial [['metacell']], top = 100)

# prioritizing metacells as the progenitor state for the metagenes 1-3
pg.mpc <- prioritize(dp, metagene = c(1, 1, 1, 0), direction = 'progenitor')

# find the differentiation paths from the progenitor toward the committed
# state of metagene 1 - 3
p2c <- differentiation.path.p2c(dp, metagene = c(1, 2, 3))

# visualizing the differentiation paths toward the committed states of metagene
# 1 - 3
dev.new(width = 7, height = 7)
plot(dp, type = 'metacell.landscape', paths = p2c$paths, 
property = pg.mpc$metacell, col.metacell = col.metacell)

}
\keyword{dpath}

