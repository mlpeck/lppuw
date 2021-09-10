
## minimum cost network flow model for phase unwrapping

#' @describeIn brcutpuw Network flow algorithm for phase unwrapping
netflowpuw <- function(phase, wts=NULL, details=FALSE, trace=1) {
  nr <- nrow(phase)
  nc <- ncol(phase)
  dx <- zernike::wrap(diff(phase))
  dy <- zernike::wrap(t(diff(t(phase))))
  d2x <- dx[,2:nc] - dx[,1:(nc-1)]
  d2y <- dy[1:(nr-1),] - dy[2:nr,]
  ndx <- sum(!is.na(dx))
  ndy <- sum(!is.na(dy))
          
  ## these weights are too simple, but this seems to work
  ## better than combining weights from each pixel that
  ## went into a given difference.
          
  if (!is.null(wts)) {
    wts <- wts/max(wts,na.rm=TRUE)
    wx <- wts[-1,]
    wy <- wts[,-1]
    wx <- wx[!is.na(dx)]
    wy <- wy[!is.na(dy)]
    obj <- c(wx, wx, wy, wy)
  } else {
    obj <- rep(1, 2*(ndx+ndy))
  }
  ## the optimal solution will cancel the map's charges
         
  charge <- -round((d2y+d2x)/(2*pi))
            
  kx <- dx
  kx[!is.na(kx)] <- 1:ndx
  ky <- dy
  ky[!is.na(ky)] <- 1:ndy
            
  ## these are the 1d indexes of the pixels that go into each 2 x 2 loop
            
  k2x1 <- kx[,1:(nc-1)][!is.na(charge)]
  k2x2 <- kx[,2:nc][!is.na(charge)]
  k2y1 <- ky[1:(nr-1),][!is.na(charge)]
  charge <- charge[!is.na(charge)]
  rm(d2x,d2y,kx,ky)
  
  nrow <- length(charge)
  ncol <- 2*(ndx+ndy)
  row <- c(-1., 1., 1., -1., 1., -1., -1., 1.)
  xvals <- rep(row, nrow)
  i_ind <- rep(1:nrow, each=8)
  j_ind <- as.numeric(rbind(k2x1, k2x2, 
                            k2x1+ndx, k2x2+ndx,
                            k2y1+2*ndx, k2y1+2*ndx+1,
                            k2y1+2*ndx+ndy, k2y1+2*ndx+ndy+1)
              )
  if (packageVersion("Matrix") >= "1.3.0") {
    mat <- Matrix::sparseMatrix(i=i_ind, j=j_ind, x=xvals, dims=c(nrow, ncol), repr="T")
  } else {
    mat <- Matrix::sparseMatrix(i=i_ind, j=j_ind, x=xvals, dims=c(nrow, ncol), giveCsparse=FALSE)
  }
  
  lpout <- rcbc::cbc_solve(obj=obj, mat=mat,
                           row_ub=charge,
                           row_lb=charge,
                           col_lb=rep(0, ncol),
                           col_ub=rep(1, ncol),
                           cbc_args = list(log=trace, verbose=trace)
                 )
  if (rcbc::solution_status(lpout) != "optimal") {
    warning("non-convergence detected")
    return(lpout)
  }
  
  dx <- rbind(dx, rep(NA, nc))
  dy <- cbind(dy, rep(NA, nr))
  soln <- lpout$column_solution
  ex <- soln[1:ndx]-soln[(ndx+1):(2*ndx)]
  ey <- soln[(2*ndx+1):(2*ndx+ndy)]-soln[(2*ndx+ndy+1):(2*ndx+2*ndy)]
  ex.m <- dx
  ey.m <- dy
  ex.m[!is.na(dx)] <- ex
  ey.m[!is.na(dy)] <- ey
  puw <- zernike::idiffpuw(phase=phase, dx=dx+2*pi*ex.m, dy=dy+2*pi*ey.m)
  if (details) {
    ex.m[ex.m==0] <- NA
    ey.m[ey.m==0] <- NA
    list(puw=puw, ex=ex.m, ey=ey.m, lpout=lpout)
  } else {
    puw
  }
                
}

