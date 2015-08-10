\name{jigsudoku}
\alias{jigsudoku}
\title{Jigsaw Sudoku solver}
\description{
  Solve a 9 x 9 jigsaw (aka squiggly) Sudoku.
}
\usage{
sudoku(givens, blocks)
}
\arguments{
  \item{givens}{Matrix of given values}
  \item{blocks}{Matrix identifying blocks}
}
\details{
  Use 0 for cells without given values.
}
\value{
  Solution as a 9 x 9 matrix.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
