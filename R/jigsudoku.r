## jigsaw sudoku solver

jigsudoku <- function(givens, block) {
	ngiven <- sum(givens > 0)
	rows <- matrix(1:9, 9, 9)
	cols <- t(rows)
	ix <- sort(as.vector(block), index.return=TRUE)$ix
	rows <- as.vector(rows)[ix]-1
	cols <- as.vector(cols)[ix]-1
	ret <- .C("jigsudoku", 
		givens = as.integer(givens), 
		ngiven = as.integer(ngiven),
		rows = as.integer(rows),
		cols = as.integer(cols),
		soln = double(9*9*9),
		retval = integer(1),
		PACKAGE = "lppuw")
	soln <- array(round(ret$soln), dim=c(9,9,9))
	if (ret$retval > 0) stop(paste("lp returned", ret$retval, sep=" "))
	nz <- which(soln>0, arr.ind=TRUE)
	sol <- matrix(0, 9, 9)
	for (i in 1:nrow(nz))
		sol[nz[i,1],nz[i,2]] <- nz[i,3]
	return(sol)
}
