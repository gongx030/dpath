library(FNN)	# knn.index
library(irlba)	# irlba
library(MASS)
library(igraph)
library(parallel)	# mclapply
library(kohonen)	# som
library(colorspace)	# hex
library(gplots)	# colorpanel
library(cluster) # pam
library(fields) # rdist

# --------------------------------------------------------------------
# [2016-07-31] Weighted NMF; When W is NULL, run a regular NMF
# --------------------------------------------------------------------
wnmf <- function(X, K, W = NULL, max.iter = 1000, threshold = 1e-3, scale.V = FALSE, verbose = TRUE){

	N <- nrow(X)
	M <- ncol(X)
	eps <- 1e-10

	if (verbose){
		cat('-------------------------------------------------------------------------------------------------\n')
		cat('Weighted non-negative matrix factorization\n')
		cat('-------------------------------------------------------------------------------------------------\n')
		cat(sprintf('number of input genes(N): %d\n', N))
		cat(sprintf('number of input cells(M): %d\n', M))
	}

	if (verbose)
		cat(sprintf('[%s] initializing metagene coefficients(V): ', Sys.time()))

	V <- irlba(X, nu = K, nv = 1)$u
	m <- kmeans(V, K)$cluster
	V <- as.matrix(as(sparseMatrix(i = m, j = 1:N, dims = c(K, N)), 'dgCMatrix') %*% X)
	V <- as.matrix(V %*% Diagonal(x = 1 / colSums(V)))

	if (verbose)
		cat('\n')

	U <- matrix(runif(N * K), nrow = N, ncol = K)	# initializing coefficients(U) randomly

	if (!is.null(W))
		XW <- X * W

	if (verbose){
		cat('-------------------------------------------------------------------------------------------------\n')
		cat(sprintf('%21.21s%5.5s%15.15s\n', '', 'iter', 'objval'))
		cat('-------------------------------------------------------------------------------------------------\n')
	}
	iter <- 1
	optval <- NULL
	while (iter <= max.iter){

		if (is.null(W)){
			U <- U * (X %*% t(V) / ( (U %*% V) %*% t(V) + eps))
			V <- V * (t(U) %*% X / (t(U) %*% (U %*% V) + eps))
		}else{
			U <- U * ((XW) %*% t(V) / ( ((U %*% V) * W) %*% t(V) + eps))
			V <- V * (t(U) %*% (XW) / (t(U) %*% ((U %*% V) * W) + eps))
		}

		if (scale.V)
			V <- V %*% diag(1 / colSums(V))

		if (iter == 1 || iter %% 10 == 0){
			Xp <- U %*% V
			if (is.null(W)){
				J <- norm(X - U %*% V, 'F')
			}else{
				J <- norm((X - U %*% V) * W, 'F')
			}
			if (verbose)
				cat(sprintf('[%s]%5.d%15.3e\n', Sys.time(), iter, J))
			optval <- rbind(optval, data.frame(iter = iter, J = J))
		}
		if (!is.null(optval) && nrow(optval) > 1){
			if (abs(optval[nrow(optval), 'J'] - optval[nrow(optval) - 1, 'J']) < threshold){
				if (verbose)
					cat(sprintf('Parameter change below threshold %.3e after %d iterations\n', threshold, iter))
				break
			}
		}

		iter <- iter + 1
	}

	list(U = U, V = V, optval = optval)

} # end of wnmf


# --------------------------------------------------------------------
# [2016-07-29] Poisson weighted non-negative matrix factorization
# --------------------------------------------------------------------
wpnmf <- function(X, K = NA, U = NULL, V = NULL, lambda0 = 0.1, max.iter = 100, threshold = 1e-3, verbose = TRUE, scale.V = TRUE){

	N <- nrow(X)	# number of genes
	M <- ncol(X)	# number of cells
	eps <- 1e-10

	if (verbose){
		cat('-------------------------------------------------------------------------------------------------\n')
		cat('Weighted Poisson non-negative matrix factorization\n')
		cat('-------------------------------------------------------------------------------------------------\n')
		cat(sprintf('number of input genes(N): %d\n', N))
		cat(sprintf('number of input cells(M): %d\n', M))
		cat(sprintf('lambda0: %.3e\n', lambda0))
	}

	if (is.null(U) || is.null(V)){
		W <- matrix(0, nrow = N, ncol = M)
		W[X > 0] <- 1
		mf <- wnmf(X, K = K, W = W, max.iter = max.iter, scale.V = scale.V, verbose = FALSE)
		U <- mf$U
		V <- mf$V
	}

	if (verbose){
		cat('-------------------------------------------------------------------------------------------------\n')
		cat(sprintf('%21.21s%5.5s%15.15s\n', '', 'iter', 'objval'))
		cat('-------------------------------------------------------------------------------------------------\n')
	}
	iter <- 1
	optval <- NULL
	P1 <- exp(ldpois(X, mu = lambda0))	# probability of being a dropout
	W <- matrix(0.5, N, M)	# initialize the weight

	while (iter <= max.iter){

		U <- U * ( ((W * X / (U %*% V + eps)) %*% t(V)) / (W %*% t(V) + eps) )
		V <- V * ( (t(U) %*% (W * X / (U %*% V + eps))) / (t(U) %*% W + eps) )
		if (scale.V)
			V <- V %*% diag(1 / colSums(V))
		P2 <- exp(ldpois(X, U %*% V))
		W <- P2 / (P1 + P2 + eps)

		if (iter == 1 || iter %% 10 == 0){
			J <- sum(W * log(P2 + eps) + (1 - W) * log(P1 + eps))
			if (verbose)
				cat(sprintf('[%s]%5.d%15.3e\n', Sys.time(), iter, J))
			optval <- rbind(optval, data.frame(iter = iter, J = J))
		}
		if (!is.null(optval) && nrow(optval) > 1){
			if (abs(optval[nrow(optval), 'J'] - optval[nrow(optval) - 1, 'J']) < threshold){
				if (verbose)
					cat(sprintf('Parameter change below threshold %.3e after %d iterations\n', threshold, iter))
				break
			}
		}
		iter <- iter + 1
	}

	list(K = K, lambda0 = lambda0, U = U, V = V, optval = optval)

} # end of wpnmf


# --------------------------------------------------------------------
# [2015-05-19] Compute metagene entropy
# --------------------------------------------------------------------
entropy2 <- function(x) -rowSums(x * log(x + .Machine$double.eps))


# --------------------------------------------------------------------
# [2016-07-30] The main function for scRNA-seq analysis 
# --------------------------------------------------------------------
dpath <- function(X, K = 5, subset.gene = NULL, subset.cell = NULL, lambda0 = 0.1, max.iter = 500, repeat.mf = 50, mc.cores = 8){

	eps <- 1e-10
	scale.V <- TRUE

	if (class(X) != 'matrix')
		stop('X must be a numeric matrix')
	
	if (any(X < 0))
		stop('there exists negative entries in expression matrix X; all entries in X must be non-negative')
	
	K <- as.integer(K)
	if (is.na(K) || K <= 0)
		stop('K must be a positive integer')
	
	if (K == 1)
		stop('K must be at least 2')

	if (K > ncol(X))
		stop('K must be smaller than the number of cells')

	if (lambda0 < 0)
		stop('lambda0 must be non-negative')

	repeat.mf <- as.integer(repeat.mf)
	if (is.na(repeat.mf) || repeat.mf <= 0)
		stop('repeat.mf must be a positive integer')

	mc.cores <- as.integer(mc.cores)
	if (is.na(mc.cores) || mc.cores<= 0)
		stop('mc.cores be a positive integer')

	N <- nrow(X)	# number of all genes
	M <- ncol(X)	# number of cells

	if (is.null(subset.gene))
		subset.gene <- rep(TRUE, N)

	if (is.null(subset.cell))
		subset.cell <- rep(TRUE, M)

	if (sum(subset.gene) <= K || sum(subset.cell) <= K)
		stop('The number of genes and cells for fitting wpnmf should be greater than K')

	cat(sprintf('-----------------------------------------------------------------------\n'))
	cat('dpath: the single cell RNA-seq analysis pipeline\n')
	cat(sprintf('-----------------------------------------------------------------------\n'))
	cat(sprintf('number of input genes=%d\n', N))
	cat(sprintf('number of input cells=%d\n', M))
	cat(sprintf('number of genes used for initializing wpnmf=%d\n', sum(subset.gene)))
	cat(sprintf('number of cells used for initializing wpnmf=%d\n', sum(subset.cell)))
	cat(sprintf('sparsity of expression matrix=%.1f%%\n', 100 * (N * M - sum(X > 0)) / (N * M)))
	cat(sprintf('number of metagenes(K)=%d\n', K))
	cat(sprintf('lambda for dropout event(lambda0)=%.3f\n', lambda0))
	cat(sprintf('number of repeated matrix factorization(repeat.mf)=%d\n', repeat.mf))
	cat(sprintf('number of cores for optimization(mc.cores)=%d\n', mc.cores))
	cat(sprintf('-----------------------------------------------------------------------\n'))

	cat(sprintf('[%s] running weighted Poisson NMF: ', Sys.time()))
	mf.list <- mclapply(1:repeat.mf, function(r){

		set.seed(r)

		# initializing U and V by a specified susbet of genes and cells
		mf <- wpnmf(log(X[subset.gene, subset.cell, drop = FALSE] + 1), K = K, lambda0 = lambda0, max.iter = max.iter, verbose = FALSE, scale.V = scale.V)

		# recovering the metagene basis for the remaining genes
		U <- matrix(NA, nrow = N, ncol = K, dimnames = list(rownames(X), NULL))
		U[subset.gene, ] <- mf$U
		if (any(!subset.gene))
			U[!subset.gene, ] <- recoverU(log(X[!subset.gene, subset.cell, drop = FALSE] + 1), mf$V, lambda0 = lambda0, max.iter = max.iter)

		# recovering the metagene coefficients for the remaining cells
		V <- matrix(NA, nrow = K, ncol = M, dimnames = list(NULL, colnames(X)))
		V[, subset.cell] <- mf$V
		if (any(!subset.cell))
			V[, !subset.cell] <- recoverV(log(X[, !subset.cell, drop = FALSE] + 1), U, lambda0 = lambda0, max.iter = max.iter, scale.V = scale.V)
		
		cat('.')

		list(U = U, V = V)
	}, mc.cores = mc.cores, mc.set.seed = FALSE)

	cat(sprintf('\n[%s] done!\n', Sys.time()))
	
	structure(list(
		X = X,	# observed expression levels
		M = M, N = N,
		K = K,	# number of metagenes							 
		repeat.mf = repeat.mf,	# number of repeated MF
		lambda0 = lambda0,
		subset.gene = subset.gene,
		subset.cell = subset.cell,
		mf.list = mf.list
	), class = 'dpath')

} # end of dpath


# --------------------------------------------------------------------
# [2015-05-21] Recover U from X and V in wpnmf
# --------------------------------------------------------------------
recoverU <- function(X, V, lambda0 = 0.1, max.iter = 100, threshold = 1e-3){

	N <- nrow(X)	# number of genes
	M <- ncol(X)	# number of cells
	K <- nrow(V)		# number of metagenes
	eps <- 1e-10

	P1 <- exp(ldpois(X, mu = lambda0))	# probability of being a dropout
	P2 <- matrix(0.5, nrow = N, ncol = M)
	W <- matrix(0.5, nrow = N, ncol = M)
	U <- matrix(runif(N * K), nrow = N, ncol = K)
	dU <- rep(0, N)	# change of U
	n <- rep(TRUE, N)
	iter <- 1

	while (iter <= max.iter){

		U0 <- U
		U[n, ] <- U[n, , drop = FALSE] * ( ((W[n, , drop = FALSE] * X[n, , drop = FALSE] / (U[n, , drop = FALSE] %*% V + eps)) %*% t(V)) / (W[n, , drop = FALSE] %*% t(V) + eps) )
		P2[n, ] <- exp(ldpois(X[n, , drop = FALSE], U[n, , drop = FALSE] %*% V))
		W[n, ] <- P2[n, , drop = FALSE] / (P1[n, , drop = FALSE] + P2[n, , drop = FALSE] + eps)
		dU[n] <- sqrt(rowSums((U[n, , drop = FALSE] - U0[n, , drop = FALSE])^2))
		n[dU < threshold] <- FALSE
		if (sum(n) == 0)
			break
		iter <- iter + 1
	}

	U
} # end of recoverU


# --------------------------------------------------------------------
# [2016-08-09] Recover V from X and U in wpnmf
# --------------------------------------------------------------------
recoverV <- function(X, U, lambda0 = 0.1, max.iter = 100, threshold = 1e-3, scale.V = TRUE){

	N <- nrow(X)	# number of genes
	M <- ncol(X)	# number of cells
	K <- ncol(U)		# number of metagenes
	eps <- 1e-10

	P1 <- exp(ldpois(X, mu = lambda0))	# probability of being a dropout
	P2 <- matrix(0.5, nrow = N, ncol = M)
	W <- matrix(0.5, nrow = N, ncol = M)
	V <- matrix(runif(K * M), nrow = K, ncol = M, dimnames = list(NULL, colnames(X)))
	dV <- rep(0, M)	# change of U
	m <- rep(TRUE, M)
	iter <- 1

	while (iter <= max.iter){

		V0 <- V
		V[, m] <- V[, m, drop = FALSE] * ( (t(U) %*% (W[, m, drop = FALSE] * X[, m, drop = FALSE] / (U %*% V[, m, drop = FALSE] + eps))) / (t(U) %*% W[, m, drop = FALSE] + eps) )

		if (scale.V)
			V[, m] <- as.matrix(V[, m, drop = FALSE] %*% Diagonal(x = 1 / colSums(V[, m, drop = FALSE])))

		P2[, m] <- exp(ldpois(X[, m, drop = FALSE], U %*% V[, m, drop = FALSE]))
		W[, m] <- P2[, m, drop = FALSE] / (P1[, m, drop = FALSE] + P2[, m, drop = FALSE] + eps)
		dV[m] <- sqrt(colSums((V[, m, drop = FALSE] - V0[, m, drop = FALSE])^2))
		m[dV < threshold] <- FALSE
		if (sum(m) == 0)
			break
		iter <- iter + 1
	}

	V	
} # end of recoverV


# --------------------------------------------------------------------
# [2016-08-07] Fitting a self-organizing map using bootstrapped cells
# --------------------------------------------------------------------
fitsom <- function(dp, k = NA, xdim = 15, ydim = 15, n.min = 10){

	som.grid <- somgrid(xdim = xdim, ydim = ydim, topo = 'hexagonal')
	M <- ncol(dp$V)	# number of cells
	K <- nrow(dp$V)
	V <- dp$V %*% diag(1 / colSums(dp$V))
	n.metacell <- som.grid$ydim * som.grid$xdim	# number of metacell
	A <- do.call('rbind', lapply(1:n.metacell, function(me) data.frame(from = which(abs(sqrt((som.grid$pts[, 'x'] - som.grid$pts[me, 'x'])^2 + (som.grid$pts[, 'y'] - som.grid$pts[me, 'y'])^2) - 1) < 1e-5),	 to = me)))	# find the neighboring metacells
	A <- sparseMatrix(i = A[, 'from'], j = A[, 'to'], dims = c(n.metacell, n.metacell))
	cat(sprintf('[%s] building %dx%d SOM:\n', Sys.time(), som.grid$xdim, som.grid$ydim))

	V2 <- cbind(as.matrix(V), diag(K))
	V2 <- V2[, sample(1:ncol(V2), n.metacell, replace = TRUE)]
	model <- som(t(V2), grid = som.grid)	# fit a SOM model by bootstraped samples
	metacell <- model$codes
	B <- sparseMatrix(i = rep(1:M, 1), j = c(knnx.index(metacell, t(V), 1)), dims = c(M, n.metacell))	# cells to metacell

	if (is.na(k)){
		# iteratively cut the SOM
		k <- 2
		while (TRUE){
			C <- sparseMatrix(i = 1:n.metacell, j = pam(metacell, k, cluster.only = TRUE), dims = c(n.metacell, k))	# metacell to clusters
			BC <- B %*% C	# cells to cluster
			Ak <- as.matrix(A); Ak[!as.matrix(C %*% t(C))] <- 0
			if (clusters(graph.adjacency(as(Ak, 'dgCMatrix')))$no > k || any(colSums(BC) < n.min))
				break
			k <- k + 1
		}
		k.max <- k - 1
	}else
		k.max <- k

	C <- sparseMatrix(i = 1:n.metacell, j = pam(metacell, k.max, cluster.only = TRUE), dims = c(n.metacell, k.max))	# metacell to clusters
	BC <- B %*% C	# cells to cluster
	center <- as.matrix(V %*% BC %*% Diagonal(x = 1 / colSums(BC)))
	dp$metacell <- metacell
	dp$cluster <- list(k = k.max, labels = max.col(BC), labels.mc = max.col(C), center = center)
	dp$cell2mc <- max.col(B)
	dp$grid <- som.grid
	dp$A <- A
	dp

} # end of fitsom


# --------------------------------------------------------------------
# [2016-07-31] Visualize the coefficients(U) and basis(V) of a list of selected genes
# --------------------------------------------------------------------
plot.dpath <- function(x, ...){

	param <- list(...)

	dp <- x
	N <- nrow(dp$X)	# number of genes
	K <- dp$K				# number of metagenes
	M <- ncol(dp$X) # number of cells
	type <- param[['type']]

	if (is.null(type)){
		cat('type must be specified as one of the following value:\n')
		cat('          markers: metagene coefficients and basis, as well as expression levels for selected marker genes\n')
		cat('     cell.cluster: cell clustering results\n')
		cat(' metagene.entropy: metagene entropy of metacells on the SOM\n')
		cat('cell.distribution: distribution of a specified group of cells on the SOM\n')
		cat('  gene.expression: gene expression level on the SOM\n')
		cat('   prioritization: top ranking genes for prioritized metacells with specified progenitor or committed state\n')
		stop('type cannot be NULL')			
	}

	# set the color for each metagene
	if (is.null(param$col.metagene))
		param$col.metagene <- rainbow(K)
	else{
		if (length(param$col.metagene) != K)
			stop('length of col.metagene must be equal to the number of metagene')
		
		is.color <- sapply(param$col.metagene, function(y) tryCatch(is.matrix(col2rgb(y)), error = function(e) FALSE))
		if (any(!is.color))
			stop('col.metagene contains non-color string')
	}

	if (type == 'markers'){

		# -----------------------------------------------------------------------------------------------
		# metagene coefficients and basis, as well as expression levels for selected marker genes
		# -----------------------------------------------------------------------------------------------

		if (is.null(param$genes)){
			gs <- rownames(dp$X)
		}else{
			is.included <- param$genes %in% rownames(dp$X) 
			if (sum(is.included) == 0)
				stop('None of supplied genes are included')
			else
				cat(sprintf('%d of %d supplied genes are included\n', sum(is.included), length(param$genes)))
			gs <- param$genes[is.included]
		}

		if (is.null(param$reorder.genes))
			param$reorder.genes <- TRUE

		D <- Reduce('+', lapply(dp$mf.list, function(mf) as.matrix(dist(t(mf$V))))) / length(dp$mf.list)
		hc <- hclust(as.dist(D))

		# panel 1: heatmap for metagenes coefficients
		x.left <- 0.05; x.right <- 0.2; y.bottom <- 0.1; y.top <- 0.9
		for (k in 1:K){
			par(fig = c(x.left + (x.right - x.left) * ((k - 1) / K), x.left + (x.right - x.left) * (k / K), y.bottom, y.top), mar = c(1, 0, 6, 0), new = ifelse(k == 1, FALSE, TRUE))
			image(t(as.matrix(dp$V[k, hc$order])), col = colorpanel(100, low = 'white', high = param$col.metagene[k]), breaks = seq(0, 1, length.out = 101), axes = FALSE)
			mtext(k, side = 3)
		}

		# panel 2: heatmap for metagene basis
		x.left <- 0.2; x.right = 1; y.bottom <- 0.025; y.top <- 0.1
		U <- as.matrix(dp$U[gs, ])
		rn <- rownames(U); U <- diag(1 / rowSums(U)) %*% U; rownames(U) <- rn

		if (param$reorder.genes){
			hc2 <- hclust(dist(U))	# clustering selected genes
			jj <- rownames(U)[hc2$order]
			U <- U[jj, ]
		}else
			jj <- rownames(U)

		for (k in 1:K){
			par(fig = c(x.left, x.right, y.bottom + (y.top - y.bottom) * ((k - 1) / K), y.bottom + (y.top - y.bottom) * (k / K)), mar = c(0, 5, 0, 5), new = TRUE)
			image(U[, k, drop = FALSE], col = colorpanel(100, low = 'white', high = param$col.metagene[k]), axes = FALSE, breaks = seq(0, 1, length.out = 101))
			mtext(k, side = 4, las = 2, line = 1)
		}

		# panel 3: barplot for median expression levels of each gene
		x.left <- 0.2; x.right <- 1; y.bottom <- 0.9; y.top <- 1
		par(fig = c(x.left, x.right, y.bottom, y.top), mar = c(0, 5, 2, 5), new = TRUE)

		X <- dp$X[jj, hc$order]	# expression matrix in log(normalized counts + 1)
		X[X < 1] <- 0
		X <- log(X + 1)
		X <- as.matrix(Diagonal(x = 1 / apply(X, 1, max)) %*% X)	# scale so that each gene's maximum expression is one
		x.median <- sapply(gs, function(g) log(mean(dp$X[g, ]) + 1))	# the mean observed expression levels

		par(lwd = 3)
		barplot(x.median[jj], names.arg = '', xaxs = 'i', ylab = 'log(count + 1)', col = 'white', cex.axis = 1.25, cex.lab = 1, main = 'Median observed count', las = 2)
		par(lwd = 1)

		# panel 4: heatmap for marker genes
		x.left <- 0.2; x.right <- 1; y.bottom <- 0.1; y.top <- 0.9
		par(fig = c(x.left, x.right, y.bottom, y.top), mar = c(1, 5, 6, 5), new = TRUE)
		image(X, col = colorpanel(100, low = 'black', high = 'green'), axes = FALSE)
		add.grid(X, col = 'black')
		col.axis <- rep('black', length(gs))
		text(y = rep(1.025, length(gs)), x = seq(0, 1, length.out = length(gs)), labels = jj, srt = 90, xpd = TRUE, col = col.axis, cex = 1.3, pos = 4, adj = c(0, 0))

	}else if (type == 'cell.cluster'){

		C <- table(dp$cluster$labels, param$cell.group)	# count how many cells belong to each cluster in each cell group
		C <- C[nrow(C):1, ]
		G <- as(sparseMatrix(i = 1:ncol(dp$V), j = dp$cluster$labels, dims = c(ncol(dp$V), dp$cluster$k)), 'dgCMatrix')
		V <- dp$V 
		V <- V %*% diag(1 / colSums(V))

		par(fig = c(0.05, 0.3, 0, 1), mar = c(3, 2, 10, 1), new = FALSE)
		image(log(t(C + 1)), col = colorpanel(100, low = 'white', high = 'lightblue'), axes = FALSE)
		axis(2, seq(0, 1, length.out = dp$cluster$k), sprintf('C%d', dp$cluster$k:1), tick = FALSE, las = 2, cex.axis = 1.5)
		axis(3, seq(0, 1, length.out = ncol(C)), colnames(C), las = 2, tick = FALSE, cex.axis = 1.5)
		coord <- expand.grid(y = seq(0, 1, length.out = dp$cluster$k), x = seq(0, 1, length.out = ncol(C)))
		mtext('# cells', 1, cex = 1.5, line = 1)
		C[C == 0] <- ''
		text(coord[, 'x'], coord[, 'y'], c(C[1:dp$cluster$k, , drop = FALSE]), cex = 1.5)

		x.left <- 0.3; x.right <- 0.7; y.bottom <- 0; y.top <- 1
		for (k in 1:K){
			par(fig = c(x.left + (x.right - x.left) * ((k - 1) / K), x.left + (x.right - x.left) * (k / K), y.bottom, y.top), mar = c(3, 0, 10, 0), new = TRUE)
			image(as.matrix(dp$cluster$center[k, dp$cluster$k:1, drop = FALSE]), col = colorpanel(100, low = 'white', high = param$col.metagene[k]), axes = FALSE, breaks = seq(0, 1, length.out = 101))
			add.grid(as.matrix(dp$cluster$center[k, dp$cluster$k:1, drop = FALSE]), col = 'black')
			if (k == round(K / 2))
				mtext('Metagene signature', 1, cex = 1.5, line = 1)
			mtext(k, 3, cex = 1.5, line = 1)
		}

		par(fig = c(0.7, 0.95, 0, 1), mar = c(3, 1, 10, 1), new = TRUE)
		barplot(sapply(split(entropy2(dp$metacell), list(dp$cluster$labels.mc)), mean)[dp$cluster$k:1], horiz = TRUE, names.arg = '', axes = FALSE, yaxs = 'i')
		axis(3, seq(0, 2, by = 0.5), cex.axis = 1.5)

	}else if (type == 'metagene.entropy'){

		# -----------------------------------------------------------------------------------------------
		# metagene entropy of metacells on the SOM
		# -----------------------------------------------------------------------------------------------

		par(mfrow = c(2, 2), mar = c(1, 2, 3, 1))
		H <- entropy2(dp$metacell)	# entropy for each metacell
		somplot(dp, property = H, topo = 'hex', main = 'Metagene entropy', cex.main = 2, col = heat.colors(100), col.border = 'white')
		z <- matrix(H, nrow = dp$grid$xdim, ncol = dp$grid$ydim)
		col <- heat.colors(100)
		colcode <- matrix(col[round( (z - min(z)) / (max(z)- min(z)) * (length(col) - 1)) + 1], nrow = nrow(z), ncol = ncol(z))
		persp(1:dp$grid$xdim, 1:dp$grid$ydim, z, theta = 10, phi = 30, xlab = '', ylab = '', zlab = 'Entropy', col = colcode, main = "Waddington's epigenetic landscape", cex.main = 2, expand = 0.25, shade = 0.25)
		somplot(dp, property = dp$metacell, topo = 'hex', main = 'Metagene distribution', cex.main = 2, col = do.call('cbind', lapply(1:K, function(k) colorpanel(100, low = 'white', high = param$col.metagene[k]))), col.border = 'white')
		par(mar = c(5, 5, 5, 2))

		if (is.null(param$cell.group))
			entropy.list <- list('All' = entropy2(t(dp$V)))
		else
			entropy.list <- split(entropy2(t(dp$V)), list(param$cell.group))
		boxplot(entropy.list, cex.axis = 2, ylab = 'Metagene entropy', cex.lab = 2)

	}else if (type == 'metacell.landscape'){

		if (is.null(param$property))
			stop('property cannot be NULL')

		if (is.null(param$col.metacell))
			param$col.metacell <- colorpanel(100, low = 'black', mid = 'white', high = 'purple')

		if (is.null(param$col.border))
			param$col.border <- 'black'

		if (is.null(param$col.cluster.border))
			param$col.cluster.border <- 'white'

		if (is.null(param$show.cluster.label))
			param$show.cluster.label <- FALSE

		somplot(dp, property = param$property, paths = param$paths, topo = 'hex', col = param$col.metacell, col.border = param$col.border, col.cluster.border = param$col.cluster.border, show.cluster.label = param$show.cluster.label, offset = 0.2)

	}else if (type == 'prioritization'){
		
		if (is.null(param$score.metacell))
			stop('score.metacell cannot be NULL')
		
		if (is.null(param$score.gene))
			stop('score.gene cannot be NULL')
		
		if (is.null(param$top))
			param$top <- 50

		if (N != length(param$score.gene))
			stop(sprintf('length of score.gene must be equal to %d', N))
		
		if (length(param$score.metacell) != nrow(dp$grid$pts))
			stop(sprintf('length of score.metacell must be equal to %d', nrow(dp$grid$pts)))
		
		gs <- names(sort(param$score.gene, decreasing = TRUE)[1:param$top])

		X <- exp(metacell.expression(dp)) - 1# observed expression levels in metacells

		# clustering the metacells by the code
		h <- hclust(dist(dp$metacell))$order

		# panel: aggregated expression levels of top ranked genes
		x.left <- 0; x.right <- 0.2; y.bottom <- 0.8; y.top <- 1
		par(fig = c(x.left, x.right, y.bottom, y.top), mar = c(1, 4, 2, 1), new = FALSE)
		somplot(dp, property = param$score.metacell, topo = 'hex', main = '', col = colorpanel(100, low = 'black', mid = 'white', high = 'purple'), col.border = 'black', show.class.label = FALSE)

		# panel: heatmap for metagenes
		x.left <- 0.05; x.right <- 0.2; y.bottom <- 0; y.top <- 0.8
		for (k in 1:K){
			par(fig = c(x.left + (x.right - x.left) * ((k - 1) / K), x.left + (x.right - x.left) * (k / K), y.bottom, y.top), mar = c(1, 0, 1, 0), new = TRUE)
			image(t(as.matrix(dp$metacell[h, k])), col = colorpanel(100, low = 'white', high = param$col.metagene[k]), axes = FALSE)
			mtext(k, side = 3)
		}

		# panel: observed expression levels at each metacell for top genes
		x.left <- 0.2; x.right <- 0.95; y.bottom <- 0; y.top <- 0.8
		par(fig = c(x.left, x.right, y.bottom, y.top), mar = c(1, 2, 1, 0), new = TRUE)
		X2 <- X[gs, h]
		X2 <- diag(1 / apply(X2, 1, max)) %*% X2
		image(X2, col = colorpanel(100, low = 'black', high = 'green'), axes = FALSE)
		add.grid(X2, col = 'black')

		# panel: metacell prioritization score
		x.left <- 0.95; x.right <- 1; y.bottom <- 0; y.top <- 0.8
		par(fig = c(x.left, x.right, y.bottom, y.top), mar = c(1, 1, 1, 1), new = TRUE)
		Z <- t(param$score.metacell[h])
		image(Z, col = colorpanel(100, low = 'black', mid = 'white', high = 'purple'), axes = FALSE)
		add.grid(Z, col = 'black')

		# panel: gene enrichment score in prioritized progenitor metacells
		x.left <- 0.2; x.right <- 0.95; y.bottom <- 0.8; y.top <- 1
		par(fig = c(x.left, x.right, y.bottom, y.top), mar = c(0, 2, 14, 0), new = TRUE)
		Z <- matrix(param$score.gene[gs], nrow = length(gs), ncol = 1)
		image(Z, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), breaks = seq(min(param$score.gene) - 1e-10, max(param$score.gene) + 1e-10, length.out = 101), axes = FALSE)
		coords <- add.grid(Z, col = 'black')

		col.axis <- rep('black', length(gs))
		text(y = rep(1.025, length(gs)), x = coords$vline[-length(coords$vline)], labels = gs, srt = 90, xpd = TRUE, col = col.axis, cex = 1, pos = 4, adj = c(0, 0))

	}else if (type == 'gene.expression'){

		if (is.null(param$genes))
			stop('genes must be specified')

		Xp <- metacell.expression(dp)	# gene expression on metacells
		for (g in param$genes){
			if (g %in% rownames(Xp)){
				somplot(dp, property = exp(Xp[g, ]) - 1, show.cluster.label = FALSE, col.cluster.border = 'white', topo = 'hex', main = g, cex.main = 2, col = colorpanel(100, low = 'black', high = 'green'), col.border = 'black')
			}
		}
	}

} # end of plot.dpath


# --------------------------------------------------------------------
# [2015-02-03] Plot SOM hex plot
# --------------------------------------------------------------------
somplot <- function(
	dp,										
	property, 
	paths = NULL,
	topo = 'hex', 
	col = heat.colors(100), 
	highlight = NULL, 
	col.highlight = 'yellow', 
	text.highlight = NULL, 
	lwd.highlight = 4,
	cex.text.highlight = 1, 
	col.text.highlight = 'black',
	col.path = 'yellow', 
	color.key = FALSE, 
	col.border = 'black', 
	col.cluster.border = 'black',
	show.cluster.label = TRUE,
	length.path.arrow = 0.075,
	lwd.path = 2, 
	offset = 0.2,
	...
){

	if (topo == 'hex')
		r <- 0.5 / cos(pi * 1/ 6)	# radius of the hex
	else
		stop(sprintf('topo %s is not supported', topo))

	n <- nrow(dp$grid$pts)	# number of metacells

	colcode <- rep('white', n)
	if (class(property) == 'numeric'){
		colcode <- col[round( (property - min(property)) / (max(property)- min(property)) * (length(col) - 1)) + 1]
	}else if (class(property) == 'matrix'){
		colcode <- do.call('cbind', lapply(1:ncol(property), function(i) col[round( (property[, i] - min(property[, i])) / (max(property[, i])- min(property[, i])) * (nrow(col) - 1)) + 1, i]))
		colcode <- sapply(1:n, function(i) {
			s <- order(property[i, ], decreasing = TRUE)[1:2]
			hex(mixcolor(0.5, hex2RGB(col2hex(colcode[i, s[1]])), hex2RGB(col2hex(colcode[i, s[2]]))))
		})
	}else
		stop('class of property should be either numeric or matrix')

	plot(0, 0, type = 'n', axes = FALSE, xlim = c(min(dp$grid$pts[, 'x']) - 0.5, max(dp$grid$pts[, 'x']) + 0.5), ylim = c(min(dp$grid$pts[, 'y']) - 0.5, max(dp$grid$pts[, 'y'] + 0.5)), xlab = '', ylab = '', ...)

	d <- do.call('rbind', lapply(1:n, function(i){	# collect the corner information for each metacell
		data.frame(
			metacell = i, 
			x = dp$grid$pts[i, 'x'] + r * cos(pi * seq(1 / 6, 2 - 1 / 6, by = 1 / 3)),
			y = dp$grid$pts[i, 'y'] + r * sin(pi * seq(1 / 6, 2 - 1 / 6, by = 1 / 3))
		)
	}))

	for (i in 1:n){
		j <- d[, 'metacell'] == i
	  polygon(d[j, 'x'], d[j, 'y'], col = colcode[i], border = col.border)
	}
	if (!is.null(dp$cluster)){

		d2 <- unique(d[, c('x', 'y')])
		nb <- summary(dp$A)
		colnames(nb) <- c('from', 'to')
		mid <- cbind(x = (dp$grid$pts[nb[, 'from'], 'x'] + dp$grid$pts[nb[, 'to'], 'x']) / 2, y = (dp$grid$pts[nb[, 'from'], 'y'] + dp$grid$pts[nb[, 'to'], 'y']) / 2)
		B <- summary(as(rdist(mid, d2) - sqrt(1/3) / 2 < 1e-10, 'dgCMatrix'))
		sp <- do.call('rbind', split(B[, 2], list(B[, 1])))
		border <- data.frame(x0 = d2[sp[, 1], 'x'], y0 = d2[sp[, 1], 'y'], x1 = d2[sp[, 2], 'x'], y1 = d2[sp[, 2], 'y'])
		i <- dp$cluster$labels.mc[nb[, 'from']] != dp$cluster$labels.mc[nb[, 'to']]
		segments(border[i, 'x0'], border[i, 'y0'], border[i, 'x1'], border[i, 'y1'], col = col.cluster.border, lwd = 2)

		if (show.cluster.label){
			sp <- split(1:n, list(dp$cluster$labels.mc))
			g <- sapply(sp, function(i) i[which.min(rowMeans(rdist(dp$grid$pts[i, , drop = FALSE])))])
			text(dp$grid$pts[g, 'x'], dp$grid$pts[g, 'y'], sprintf('C%d', 1:max(dp$cluster$labels.mc)), cex = 1.5)
		}
		
	}

	if (!is.null(highlight)){
		for (i in highlight){
		  polygon(
				dp$grid$pts[i, 'x'] + r * cos(pi * seq(1 / 6, 2 + 1 / 6, by = 1 / 3)),
				dp$grid$pts[i, 'y'] + r * sin(pi * seq(1 / 6, 2 + 1 / 6, by = 1 / 3)),
				col = colcode[i], border = col.highlight, lwd = lwd.highlight
			)
			if (!is.null(text.highlight))
				text(dp$grid$pts[highlight, 'x'], dp$grid$pts[highlight, 'y'], text.highlight, cex = cex.text.highlight, col = col.text.highlight)
		}
	}

	if (!is.null(paths)){
		z <- do.call('rbind', lapply(paths, function(p) data.frame(from = p[1:(length(p) - 1)], to = p[2:length(p)])))
		from.x <- dp$grid$pts[z[, 'from'], 'x']
		from.y <- dp$grid$pts[z[, 'from'], 'y']
		to.x <- dp$grid$pts[z[, 'to'], 'x']
		to.y <- dp$grid$pts[z[, 'to'], 'y']
		arrows(from.x + (to.x - from.x) * offset, from.y + (to.y - from.y) * offset, to.x + (from.x - to.x) * offset, to.y + (from.y - to.y) * offset, length = length.path.arrow, lwd = lwd.path, col = col.path)
	}

} # end of somplot


# --------------------------------------------------------------------
# [2016-05-22] Reorder metagenes from repeated runs
# --------------------------------------------------------------------
reorder.dpath <- function(x, ...){

	dp <- x
	K <- dp$K
	repeat.mf <- length(dp$mf.list)
	cat(sprintf('[%s] reordering and merging metagene basis/coefficients\n', Sys.time()))
	V.list <- lapply(dp$mf.list, function(mf) mf$V)	# a list of metagene coefficient of repeated runs

	# for each pair of repeated run, find the most similar metagenes
	C <- matrix(NA, nrow = repeat.mf, ncol = repeat.mf)
	for (i in 1:repeat.mf){
		for (j in 1:repeat.mf){
			g <- max.col(cor(as.matrix(t(V.list[[i]])), as.matrix(t(V.list[[j]]))))
			C[i, j] <- length(unique(g))
		}
	}

	best <- which.max(rowSums(C == K))	# the repeated run that has most distinctive metagene
	reordered.mf <- rep(FALSE, repeat.mf)
	for (i in 1:repeat.mf){
		g <- max.col(cor(as.matrix(t(V.list[[best]])), as.matrix(t(V.list[[i]]))))
		if (length(unique(g)) == K){
			V.list[[i]] <- V.list[[i]][g, ]
			dp$mf.list[[i]]$U <- dp$mf.list[[i]]$U[, g]
			reordered.mf[i] <- TRUE
		}
	}

	# Compute the mean U and V of repeated runs
	U <- Reduce('+', lapply(dp$mf.list[C[best, ] == K], function(mf) mf$U)) / sum(C[best, ] == K)	# combined coefficients(U)
	V <- Reduce('+', V.list[C[best, ] == K]) / sum(C[best, ] == K)	# combined basis(V)
	colnames(V) <- colnames(dp$X)
	dp$U <- U; dp$V <- V; dp$reordered.mf <- reordered.mf
	dp

} # end of reorder.dpath


# --------------------------------------------------------------------
# [2015-05-27] Prioritize metacells
# --------------------------------------------------------------------
prioritize <- function(dp, direction = NULL, metagene = NULL, beta = 0.9){

	if (is.null(metagene))
		stop('metagene cannot be NULL')
	
	if (!direction %in% c('progenitor', 'committed'))
		stop('direction must be either progenitor or committed')
	
	if (beta <= 0 || beta >= 1)
		stop('beta must be greater than 0 and smaller than 1')
	
	K <- dp$K	# number of metagenes
	H <- dp$grid$xdim * dp$grid$ydim	# number of metacells
	N <- nrow(dp$U)	# number of genes

	if (!is.vector(metagene) || length(metagene) != K || any(metagene < 0))
		stop(sprintf('metagene must be a %d-length non-negative vector', K))

	metagene <- metagene / sum(metagene)

	cat(sprintf('[%s] prioritizing %s metacells for metagene signature:\n%s\n', Sys.time(), direction, paste(sprintf('MG%d:%.3f', 1:K, metagene), collapse = '; ')))
	A <- transit.prob(dp$A, dp$metacell, direction = direction)	# transition probability between metacells
	u <- rwr.mcmg(A, dp$metacell, metagene, beta)$metacell

	# scale so that for each bootstrap sample the sum of metacell rank equal to zero
	rk <- u - mean(u)

	cat(sprintf('[%s] computing the gene enrichment score:\n', Sys.time()))
	W <- metacell.expression(dp)	# observed metacell expression
	W <- W + 1	# adding a pseudo count to the observed metacell expression
	w <- 1 / rowSums(W)
	w[is.infinite(w)] <- 0	
	W <- Diagonal(x = w) %*% W	# scale the sum of observed expression in all metacell to one
	score <- (W %*% rk)[, 1]
	list(gene = score, metacell = u)

} # end of prioritize


# --------------------------------------------------------------------
# [2015-05-19] Preparing the transition probability matrix between metacells
# --------------------------------------------------------------------
transit.prob <- function(A, V, direction = 'committed', p0 = 1e-3){
	K <- ncol(V)
	V <- V + .Machine$double.eps
	a <- summary(A)
	h <- entropy2(V)	# entropy for each metacell
	a[, 3] <- rowSums((V[a[, 1], ] - V[a[, 2], ])^2)
	a[, 3] <- exp(-a[, 3] / 2)	# converting distance to similarity
	a[h[a[, 1]] < h[a[, 2]], 3] <- p0	# set a low probability to reverse entropy links
	A <- sparseMatrix(i = a[, 1], j = a[, 2], x = a[, 3], dims = c(nrow(V), nrow(V)))
	if (direction == 'progenitor')
		A <- t(A)
	w <- 1 / rowSums(A); w[is.infinite(w)] <- 0
	A <- Diagonal(x = w) %*% A
	A
} # end of transit.prob


# --------------------------------------------------------------------
# [2015-04-30] Random walk with restart on a metagene-metacell graph
# --------------------------------------------------------------------
rwr.mcmg <- function(A, V, metagene = NULL, beta = 0.9){
	
	K <- ncol(V)
	M <- nrow(V)
	
	G.mcmg <- V	# metacell -> metagene graph

	# metagene -> metacell graph
	G.mgmc <- t(V)
	w <- 1 / rowSums(G.mgmc)
	w[is.infinite(w)] <- 0
	G.mgmc <- Diagonal(x = w) %*% G.mgmc

	G.mg <- Diagonal(n = K)

	w <- 1/ rowSums(A); w[is.infinite(w)] <- 0; A <- Diagonal(x = w) %*% A

	G <- rBind(
		cBind(beta * G.mg, (1 - beta) * G.mgmc), 
		cBind((1 - beta) * G.mcmg, beta * A)
	)
	v <- rep(0, K + M)
	v[1:K] <- metagene	# setting the start metagene
	u <- page.rank(graph.adjacency(as(G, 'dgCMatrix'), mode = 'directed', weighted = TRUE), personalized = v)$vector
	u.mg <- u[1:K]
	u.mc <- u[(K + 1):(K + M)]	# extract the probability on metacell
	list(metagene = u.mg / sum(u.mg), metacell = u.mc / sum(u.mc))

} # end rwr.mcmg


# --------------------------------------------------------------------
# [2014-12-05] Get observed expression in metacell
# --------------------------------------------------------------------
metacell.expression <- function(dp, k = 10, min.tpm = 0){
	nn <- knnx.index(t(dp$V), dp$metacell, k = k)
	A <- sparseMatrix(i = rep(1:nrow(dp$metacell), k), j = c(nn), x = rep(1 / k, length(nn)), dims = c(nrow(dp$metacell), ncol(dp$V)))
	X <- as.matrix(log(dp$X + 1) %*% t(A))
	X[X <= log(1 + min.tpm)] <- 0
	X
} # end of meteacell.expression


# --------------------------------------------------------------------
# [2015-04-21] Predict the shortest differentiation pathways from predicted 
# progenitor states to committed states
# --------------------------------------------------------------------
differentiation.path.p2c <- function(dp, metagene = NULL){
	
	if (is.null(metagene))
		stop('metagene must be specified')
	
	K <- nrow(dp$V)	# number of metagenes
	H <- nrow(dp$grid$pts)	# number of metacells

	if (any(metagene < 1) || any(metagene > K))
		stop(sprintf('metagene must be between 1 and %d', K))

	# prioritize the metacells for committed states
	committed <- lapply(metagene, function(k){
		v <- rep(0, K); v[k] <- 1
		prioritize(dp, metagene = v, direction = 'committed')
	})

	# prioritze the metacells for common progenitor states
	v <- rep(0, K); v[metagene] <- 1; v <- v / sum(v)
	progenitor <- prioritize(dp, metagene = v, direction = 'progenitor')

	mcc <- lapply(committed, function(x) which.max(x$metacell))
	mcp <- which.max(progenitor$metacell)

	A <- transit.prob(dp$A, dp$metacell, direction = 'progenitor')	# transition probability between metacells
	A <- summary(A)
	A[, 3] <- -log(A[, 3])	# convert similarty to distance
	A <- sparseMatrix(i = A[, 1], j = A[, 2], x = A[, 3], dims = c(H, H))
	A <- graph.adjacency(A, mode = 'directed', weighted = TRUE)
	paths <- lapply(mcc, function(to) as.vector(get.all.shortest.paths(A, from = mcp, to = to)$res[[1]]))
	
	list(mcp = mcp, paths = paths, A = A)

}


# --------------------------------------------------------------------
# [2014-06-27] log density of Poission distribution
# --------------------------------------------------------------------
ldpois <- function(k, mu){

	mu <- mu + .Machine$double.eps

	if (any(k <= -1))
		stop('k must be greater than -1')
	
	k * log(mu) - mu - lgamma(k + 1)

} # end of ldpois


# --------------------------------------------------------------------
# [2016-12-02] Add gridlines to a heatmap
# --------------------------------------------------------------------
add.grid <- function(x, col = 'gray', lwd = 0.1, v = TRUE, h = TRUE){

	if (!is.matrix(x))
		stop('x must be a matrix')

	if (ncol(x) == 1)
		ch <- 1
	else
		ch <- 1 / (ncol(x) - 1)

	if (nrow(x) > 1)
		cw <- 1 / (nrow(x) - 1)
	
	vline <- NULL
	hline <- NULL

	if (v){
		if (nrow(x) == 1)
			vline <- c(-1, 1)
		else
			vline <- seq(-cw / 2, 1 + cw / 2, length.out = nrow(x) + 1)
		abline(v = vline, col = col, lwd = lwd, xpd = FALSE)
	}

	if (h){
		if (ncol(x) == 1)
			hline <- c(-1, 1)
		else
			hline <- seq(-ch / 2, 1 + ch / 2, length.out = ncol(x) + 1)
		abline(h = hline, col = col, lwd= lwd, xpd = FALSE)	
	}
	list(vline = vline, hline = hline)	
}
