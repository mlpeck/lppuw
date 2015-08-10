\name{netflowpuw}
\alias{netflowpuw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Minimum cost network flow model for two dimensional phase unwrapping.}
\description{
  Unwraps a 2D phase map using a minimum cost network flow model
  originally proposed by M. Constantini
  (\url{http://www.geo.unizh.ch/rsl/fringe96/papers/costantini/}).
}
\usage{
netflowpuw(phase, wts=NULL, details = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phase}{Matrix of wrapped phase values.}
  \item{wts}{Matrix of quality values for weighting.}
  \item{details}{Boolean: return details of solution?}
}
\value{
  
  If \code{details} a list with items
  \item{puw}{Unwrapped phase}
  \item{ex}{Phase discontinuities - x}
  \item{ey}{Phase discontinuities - y}
  \item{objval}{Objective function value}
  
  otherwise the unwrapped phase
}
\references{\url{http://www.geo.unizh.ch/rsl/fringe96/papers/costantini/}}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{The current version minimizes the total absolute value
  of the discontinuities.
  
}
\seealso{\code{\link{brcutpuw}}}
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
  ph.uw <- netflowpuw(ph)
  
  screen(2)
  
  # contour plot of unwrapped phase
  
  plot(ph.uw, axes=FALSE)
  close.screen(all.screens=T)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
