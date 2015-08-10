sudoku <- function(givens) {
	ngiven <- sum(givens > 0)
	ret <- .C("sudoku", 
		givens = as.integer(givens), 
		ngiven = as.integer(ngiven),
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
