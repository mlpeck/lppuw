\name{brcutpuw}
\alias{brcutpuw}
\title{BRanch CUT algorithm for phase unwrapping}
\description{
  Branch cut algorithm for two dimensional phase unwrapping.
}
\usage{
  brcutpuw(phase, pen=0, details=FALSE)
}
\arguments{
  \item{phase}{A matrix of wrapped phase values}
  \item{pen}{Penalty to be added to ``cost'' of connecting a residue to an edge.}
  \item{details}{Return a map of the branch cuts?}
}
\value{
  if \code{details==TRUE}, a list consisting of:
  \item{puw}{A matrix of class \code{"\link{pupil}"} with the unwrapped wavefront.}
  \item{bcuts}{map of the branch cuts}
  
  If \code{!details} the unwrapped wavefront.
  
}
\references{
  Ghiglia, D.C., and Pritt, M.D., 1998, \bold{Two-Dimensional Phase Unwrapping},
  New York: Wiley & Sons, Inc., ISBN 0-471-24935-1.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{\code{\link{netflowpuw}}}
\note{
  The current version now solves a variation of the assignment problem
  to minimize the length of branch cuts. If \code{pen == 0} the total
  length of all branch cuts is minimized subject to the constraint
  that each residue connects to exactly one residue of the opposite
  charge, or to the closest edge point.
}
\section{Warning}{
  It is assumed that at least one row and one column on all sides of the
  input \code{phase} matrix is masked, i.e. has values \code{NA}.
  No check is performed for this and the function will return an error
  if the condition is not met.
}
\examples{
  require(zernike)
  
  # defined area must be surrounded by at least one
  # row/column of NA's
  
  ph <- matrix(NA, 27, 27)
  ph[2:26,2:26] <- runif(25^2, -pi, pi)
  
  x11(width=12, height=6)
  par(mar=c(0,0,1,0))
  split.screen(figs=c(1,2))
  screen(1)
  
  # show the wrapped phase map with corrupted pixels
  
  mtext(rmap(ph, plot=TRUE))
  ph.uw <- brcutpuw(ph)
  
  screen(2)
  
  # contour plot of unwrapped phase
  
  plot(ph.uw, axes=FALSE)
  close.screen(all.screens=T)
}
\keyword{mathematics}

