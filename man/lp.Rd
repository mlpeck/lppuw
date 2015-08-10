\name{lp}
\alias{lp}
\title{Solve a linear program}
\description{
  Solves a generic linear program with continuous
  (not integer) variables using the lp\_solve API
  (see references).
}
\usage{
lp(dir = "min", objective, constr.mat, constr.dir, rhs)
}
\arguments{
  \item{dir}{direction of optimization: "min" or "max"}
  \item{objective}{Vector of objective function coefficients}
  \item{constr.mat}{Matrix of constraint coefficients}
  \item{constr.dir}{Vector of constraint directions}
  \item{rhs}{Vector of right hand side values}
}
\details{
  \code{constr.dir} is a vector of strings. Recognized values are
  "<", "<=", "=", "==", ">", ">=".
  Anything else generates an error.
}
\value{
  \item{objval}{Objective function value}
  \item{soln}{Solution vector}
}
\references{\url{https://sourceforge.net/projects/lpsolve/}}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{At present no check is made that lp\_solve actually returned
  an optimal solution.
}
\examples{
## A toy lp

lp("max", c(143,160), matrix(c(120,110,1,210,30,1),3,2),
  rep("<", 3), c(15000,4000,75))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
