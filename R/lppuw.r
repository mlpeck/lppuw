## Package lppuw
## Routines for linear program based phase unwrapping.
## Also a simple generic lp routine
## Incorporates code for lp_solve (http://lpsolve.sourceforge.net/5.5/).

## Author: M.L. Peck (mpeck1@ix.netcom.com)
## Language: R (http://www.r-project.org/)
## Copyright (c) 2008, M.L. Peck
## Last mod: July 2008


## This is free software and is offered with no warranty whatsoever.
## License is GNU GPL (http://www.gnu.org/licenses/gpl.html).

##.First.lib <- function(libname, pkgname)
##{
##    library.dynam("lppuw", pkgname, libname)
##}

## generic linear program

lp <- function(dir = "min", objective, constr.mat, constr.dir, rhs) {
	if (dir == "min")
		dir <- 0
	else
		dir <- 1
	nvars <- length(objective)
	nconstr <- nrow(constr.mat)
	objective = c(0, objective)
	
	constrtype <- rep(-1, length(constr.dir))
	constrtype[constr.dir == "<" | constr.dir == "<="] <- 1
	constrtype[constr.dir == "=" | constr.dir == "=="] <- 3
	constrtype[constr.dir == ">" | constr.dir == ">="] <- 2
	if(any(constrtype == -1))
		stop("Unknown constraint direction found\n")
	
	ret <- .C("gen_lp",
	   dir = as.integer(dir),
	   nvars = as.integer(nvars),
	   objective = as.double(objective),
	   nconstr = as.integer(nconstr),
	   constr = as.double(constr.mat),
	   constrtype = as.integer(constrtype),
	   rhs = as.double(rhs),
	   objval = double(1),
	   soln = double(nvars),
	   PACKAGE = "lppuw")
	return(list(objval=ret$objval, soln=ret$soln))
}

## minimum cost network flow model for phase unwrapping
	
netflowpuw <- function(phase, wts=NULL, details=FALSE) {
	require(zernike)
	cat("Warning: this routine can take a very long time\n")
	flush.console()
	
	nr <- nrow(phase)
	nc <- ncol(phase)
	dx <- wrap(diff(phase))
	dy <- wrap(t(diff(t(phase))))
	d2x <- dx[,2:nc] - dx[,1:(nc-1)]
	d2y <- dy[1:(nr-1),] - dy[2:nr,]
	
	## these weights are too simple, but this seems to work
	## better than combining weights from each pixel that
	## went into a given difference.
	
        if (!is.null(wts)) {
            wts <- wts/max(wts,na.rm=TRUE)
            wx <- wts[-1,]
            wy <- wts[,-1]
            wx <- wx[!is.na(dx)]
            wy <- wy[!is.na(dy)]
            objective = c(0, wx, wx, wy, wy)
        }
        else objective = NULL	
	
	## the optimal solution will cancel the map's charges
	
	charge <- -round((d2y+d2x)/(2*pi))
	
	ndx <- sum(!is.na(dx))
	ndy <- sum(!is.na(dy))
	kx <- dx
	kx[!is.na(kx)] <- 1:ndx
	ky <- dy
	ky[!is.na(ky)] <- 1:ndy
	
	## these are the 1d indexes of the pixels that go into each 2 x 2 loop
	
	k2x1 <- kx[,1:(nc-1)][!is.na(charge)]
	k2x2 <- kx[,2:nc][!is.na(charge)]
	k2y1 <- ky[1:(nr-1),][!is.na(charge)]
	nloops <- length(k2x1)
	rm(d2x,d2y,kx,ky)
	
	lpout <- .C("netflow", 
		ndx = as.integer(ndx), 
		ndy = as.integer(ndy),
		nconstr = as.integer(nloops),
		k2x1 = as.integer(k2x1), 
		k2x2 = as.integer(k2x2),
		k2y1 = as.integer(k2y1), 
		charge = as.integer(charge[!is.na(charge)]),
		objective = as.double(objective), 
		objval = double(1), 
		soln = double(2*ndx+2*ndy),
		ret = integer(1),
		PACKAGE = "lppuw")
	if (lpout$ret != 0) stop(paste("lp_solve API returned", lpout$ret, sep=" "))
	lpout <- list(objval = lpout$objval, soln = lpout$soln)
	dx <- rbind(dx, rep(NA, nc))
	dy <- cbind(dy, rep(NA, nr))
	ex <- lpout$soln[1:ndx]-lpout$soln[(ndx+1):(2*ndx)]
	ey <- lpout$soln[(2*ndx+1):(2*ndx+ndy)]-lpout$soln[(2*ndx+ndy+1):(2*ndx+2*ndy)]
	ex.m <- dx
	ey.m <- dy
	ex.m[!is.na(dx)] <- ex
	ey.m[!is.na(dy)] <- ey
	puw <- idiffpuw(phase=phase, dx=dx+2*pi*ex.m, dy=dy+2*pi*ey.m)
	if (details) {
		ex.m[ex.m==0] <- NA
		ey.m[ey.m==0] <- NA
		return(list(puw=puw, ex=ex.m, ey=ey.m, objval=lpout$objval))
	}
	else return(puw)
	
}

##
##  Branch cut algorithm.
## This now solves a variant of the assignment problem to minimize
## the total length of all branch cuts.
##
## note: parameter pen is a penalty added to distances from charges to edges.
##		 I'm not sure this is really useful. Seems to work fine with the default
##		 of 0.
##

	
brcutpuw <- function(phase, pen=0, details=FALSE) {

	require(zernike)

	## distance between points specified by their x,y coordinates.
	
	dmat <- function(p1, p2) {
		fn <- function(x,y) abs(y-x)
		xleg <- outer(p1[,1], p2[,1], FUN=fn)
		yleg <- outer(p1[,2], p2[,2], FUN=fn)
		return(xleg+yleg)
	}
	
	## Make a branch cut
	
	brcut <- function(startp, endp, mask) {
		xl <- endp[1]-startp[1]
		yl <- endp[2]-startp[2]
		if (abs(yl) <= abs(xl)) {
			slope <- yl/xl
			ix <- startp[1]:endp[1]
			iy <- startp[2] + round(slope*(ix-startp[1]))
		} else {
			slope <- xl/yl
			iy <- startp[2]:endp[2]
			ix <- startp[1] + round(slope*(iy-startp[2]))
		}
		mask[cbind(ix,iy)] <- NA
		return(mask)
	}
	
	res <- rmap(phase)
	## no residues, so a single call to idiffpuw is all we need
	if (sum(abs(res),na.rm=TRUE) == 0) 
		return(idiffpuw(phase,ucall=TRUE))
	

	nr <- nrow(phase)
	nc <- ncol(phase)
	dx <- wrap(.fdiff(phase))
	dy <- wrap(t(.fdiff(t(phase))))
	
	cutlen <- NULL

	mask <- matrix(1,nr,nc)
	mask[is.na(phase)] <- NA
	
	## positive charges
	chp <- which(res==1, arr.ind=TRUE)
	ncp <- nrow(chp)
	## negative charges
	chm <- which(res== -1, arr.ind=TRUE)
	ncm <- nrow(chm)
	## get the boundary points - ones one pixel outside the area of defined phase.
	## step in
	ptsb <- which(is.na(phase) & (!is.na(rbind(phase[-1,], rep(NA, nc)))), arr.ind=TRUE)
	ptsb <- rbind(ptsb, which(is.na(phase) & (!is.na(cbind(phase[,-1], rep(NA, nr)))), arr.ind=TRUE))
	## step out
	ptsc <- which((!is.na(phase)) & is.na(rbind(phase[-1,], rep(NA, nc))), arr.ind=TRUE)
	ptsc[,1] <- ptsc[,1]+1
	ptsb <- rbind(ptsb, ptsc)
	ptsc <- which((!is.na(phase)) & is.na(cbind(phase[,-1], rep(NA, nr))), arr.ind=TRUE)
	ptsc[,2] <- ptsc[,2]+1
	ptsb <- rbind(ptsb, ptsc)
	
	
	if ((ncp>0) & (ncm>0)) {
	
	## distances between residues
	dpm <- dmat(chp, chm)
	
	## distance from each positive charge to nearest edge
	dce <- dmat(chp, ptsb)
	ipb <- apply(dce, 1, which.min)
	dpe <- apply(dce, 1, min)+pen
	pe <- matrix(ptsb[ipb,], ncp, 2)
	
	## distance from each negative charge to nearest edge
	dce <- dmat(chm, ptsb)
	ipb <- apply(dce, 1, which.min)
	dme <- apply(dce, 1, min)+pen
	me <- matrix(ptsb[ipb,], ncm, 2)
	
	## set up the LP. This is almost the assignment problem.
	
	obj <- c(0, as.vector(dpm), dpe, dme)
	
	ret <- .C("mod_assign", 
		ncp = as.integer(ncp),
		ncm = as.integer(ncm), 
		objective = as.double(obj),
		objval = double(1), 
		soln = double(length(obj)-1),
		retval = integer(1),
		PACKAGE = "lppuw")
	if (ret$retval != 0) stop(paste("lp_solve API returned", ret$retval, sep=" "))
	lpsol <- ret$soln
	cutlen <- ret$objval
	dpcuts <- matrix(round(lpsol[1:(ncp*ncm)]), ncp, ncm)
	pecuts <- round(lpsol[(ncp*ncm+1):(ncp*ncm+ncp)])
	mecuts <- round(lpsol[(ncp*ncm+ncp+1):(ncp*ncm+ncp+ncm)])
	dpcuts <- which(dpcuts==1, arr.ind=TRUE)
	i <- 1
	while (i <= nrow(dpcuts)) {
		startp <- chp[dpcuts[i,1],]
		endp <- chm[dpcuts[i,2],]
		mask <- brcut(startp, endp, mask)
		i <- i+1
	}
	pecuts <- which(pecuts==1)
	i <- 1
	while (i <= length(pecuts)) {
		startp <- chp[pecuts[i],]
		endp <- pe[pecuts[i],]
		mask <- brcut(startp, endp, mask)
		i <- i+1
	}
	mecuts <- which(mecuts==1)
	i <- 1
	while (i <= length(mecuts)) {
		startp <- chm[mecuts[i],]
		endp <- me[mecuts[i],]
		mask <- brcut(startp, endp, mask)
		i <- i+1
	}
	}
	## if there are only plus or minus charges connect to nearest edge (should be rare case, but...)
	else if (ncp > 0) {
	## distance from each positive charge to nearest edge
	dce <- dmat(chp, ptsb)
	ipb <- apply(dce, 1, which.min)
	dpe <- apply(dce, 1, min)+pen
	pe <- matrix(ptsb[ipb,],ncp,2)
	for (i in 1:ncp) {
		startp <- chp[i,]
		endp <- pe[i,]
		mask <- brcut(startp, endp, mask)
	}
	}
	else if (ncm > 0) {	
	## distance from each negative charge to nearest edge
	dce <- dmat(chm, ptsb)
	ipb <- apply(dce, 1, which.min)
	dme <- apply(dce, 1, min)+pen
	me <- matrix(ptsb[ipb,],ncm,2)
	for (i in 1:ncm) {
		startp <- chm[i,]
		endp <- me[i,]
		mask <- brcut(startp, endp, mask)
	}
	}
	
	## call idiffpuw to unwrap all unmasked pixels
	uwum <- idiffpuw(phase, mask, ucall=FALSE, dx=dx, dy=dy)
	puw <- uwum$puw
	uw <- uwum$uw
	uw[is.na(phase)] <- NA
	## fill in whatever was left wrapped
	## This inverts the logic of fill algorithm in idiffpuw
	## todo initially contains a list of pixels that haven't been unwrapped
	## that should be. We look for neighbors that have been unwrapped
	## and unwrap from them, removing the newly unwrapped pixels from
	## the todo list.
	## This should work because branch cuts are just a pixel wide and should
	## adjoin unmasked regions. It's possible that cuts could completely
	## enclose an area, in which case it will still be filled
	## with most likely erroneous values as we fill through the cuts.
	
	todo <- which(!uw)
	kn <- c(-1, 1, -nr, nr)
	while (length(todo) > 0) {
		for (i in 1:4) {
			b <- todo + kn[i]
			subs <- which(uw[b])
			puw[todo[subs]] <- puw[b[subs]] + switch(i, dx[b[subs]], -dx[todo[subs]],
			 dy[b[subs]], -dy[todo[subs]])
			uw[todo[subs]] <- TRUE
			todo <- todo[-subs]
			if (length(todo)==0) break
		}
	}
	if (details) {
		bcuts <- matrix(NA, nr, nc)
		bcuts[is.na(mask) & (!is.na(phase))] <- 1
		return(list(puw=puw/(2*pi), cutlen=cutlen, bcuts=bcuts))
	} else
		return(puw/(2*pi))
}
