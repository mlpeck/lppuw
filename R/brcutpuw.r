##
##  Branch cut algorithm.
## This solves a variant of the assignment problem to minimize
## the total length of all branch cuts.
##
## note: parameter pen is a penalty added to distances from charges to edges.
##		 I'm not sure this is really useful. Seems to work fine with the default
##		 of 0.
##

#' Linear program based phase unwrapping
#'
#' `brcutpuw` implements a branch cut algorithm for phase unwrapping
#'
#' @details This package implements two distinct algorithms for two dimensional
#'   phase unwrapping that can be set up and solved as general linear programs.
#'   The linear programs are solved using the LP solver `CBC` from [COIN-OR](https://github.com/coin-or/Cbc)
#'   using the interface provided by the package [rcbc](https://dirkschumacher.github.io/rcbc/index.html).
#'
#' @param phase matrix with phase map to be unwrapped (phases are in radians)
#' @param pen penalty for making a branch cut from a residue to an edge (`brcutpuw` only)
#' @param wts matrix of weights for cost function with same dimension as `phase` (`netflowpuw` only)
#' @param details Return some details of the solution?
#' @param trace Send some info from the LP solver to the console?
#'
#' @return The unwrapped wavefront in units of fringes.
#'   If `details` is `TRUE` additional details of the solution
#'   are returned in a named list with the first member
#'   `puw` containing the unwrapped wavefront.
#'
#' @seealso There is a function with the same name `brcutpuw` in package [zernike].
#'
#' @describeIn brcutpuw Branch cut algorithm for phase unwrapping
#'
#' @section Note:
#'   According to the documentation for `rcbc` different levels of detail
#'   from the LP solver can be printed with trace levels up to 15,
#'   however the same output seems to be returned for all values.
#'   Setting `trace=0` will produce silent output, which may not
#'   be advisable since these can take some time to run.
#'
#'   The value of trace doesn't matter when running in a Windows "GUI"
#'   console window or Rstudio because the CBC output isn't passed. If
#'   you want to see some output run in a terminal window instead.
#'
#' @examples
#'   data("phasemaps", package="lppuw")
#'   mtext(zernike::rmap(phi, plot=TRUE))
#'   wf.bc <- brcutpuw(phi)
#'   wf.nf <- netflowpuw(phi, mod)
#'   X11()
#'   zernike::plot.pupil(wf.nf, col=zernike::rygcb(400))
#'   cat("Summary of the difference between the two unwrapped wavefronts:\n")
#'   zernike::summary.pupil(wf.nf - wf.bc)
#'
brcutpuw <- function(phase, pen=0, details=FALSE, trace=1) {
          
  ## distance between points specified by their x,y coordinates.
          
  dmat <- function(p1, p2) {
    fn <- function(x,y) abs(y-x)
    xleg <- outer(p1[,1], p2[,1], FUN=fn)
    yleg <- outer(p1[,2], p2[,2], FUN=fn)
    xleg+yleg
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
    mask
  }
  
  res <- zernike::rmap(phase)
  ## no residues, so a single call to idiffpuw is all we need
  if (sum(abs(res),na.rm=TRUE) == 0) {
    return(zernike::idiffpuw(phase,ucall=TRUE))
  }
            
    
  nr <- nrow(phase)
  nc <- ncol(phase)
  dx <- zernike::wrap(zernike::.fdiff(phase))
  dy <- zernike::wrap(t(zernike::.fdiff(t(phase))))
    
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
      
    obj <- c(as.vector(dpm), dpe, dme)
    nrow <- ncp + ncm
    ncol <- ncp*ncm + ncp + ncm
    row_lb <- rep(1, nrow)
    i_ind <- c(rep(1:ncp, each=ncm+1), ncp+rep(1:ncm, each=ncp+1))
    j_mp <- t(matrix(1:ncp, ncp, ncm)) + matrix(ncp*(0:(ncm-1)), ncm, ncp)
    j_p <- rbind(j_mp, (ncm*ncp)+(1:ncp))
    j_p <- as.vector(j_p)
    j_m <- rbind(t(j_mp), (ncm*ncp + ncp)+(1:ncm))
    j_m <- as.vector(j_m)
    j_ind <- c(j_p, j_m)
    mat <- Matrix::sparseMatrix(i=i_ind, j=j_ind, dims=c(nrow, ncol), giveCsparse=FALSE)
    lpout <- rcbc::cbc_solve(obj=obj, mat=mat,
                            row_lb = row_lb,
                            row_ub = row_lb,
                            col_lb = rep(0, ncol),
                            col_ub = rep(1, ncol),
                            cbc_args = list(log=trace, verbose=trace)
                          )
    if (rcbc::solution_status(lpout) != "optimal") {
      warning("non-convergence detected")
      return(lpout)
    }
      
    lpsol <- lpout$column_solution
    cutlen <- lpout$objective_value
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
  } else if (ncp > 0) {
      ## if there are only plus or minus charges connect to nearest edge (should be rare case, but...)

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
  } else if (ncm > 0) {
    
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
  
  uwum <- zernike::idiffpuw(phase, mask, ucall=FALSE, dx=dx, dy=dy)
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
    list(puw=puw/(2*pi), cutlen=cutlen, bcuts=bcuts, 
         obj=obj, mat=mat, lpout=lpout)
  } else {
    puw/(2*pi)
  }
}
