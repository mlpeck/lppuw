\name{sudoku}
\alias{sudoku}
\title{Sudoku solver}
\description{
  Solve a standard 9 x 9 Sudoku.
}
\usage{
sudoku(givens)
}
\arguments{
  \item{givens}{Matrix of given values}
}
\details{
  Use 0 for cells without given values.
}
\value{
  Solution as a 9 x 9 matrix.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\examples{
puzzle <- matrix(c(1, 0, 0, 0, 0, 0, 8, 0, 6,
 0, 0, 8, 0, 0, 6, 0, 0, 0, 
 0, 4, 0, 7, 8, 0, 0, 5, 0, 
 4, 0, 0, 0, 6, 0, 0, 0, 8, 
 7, 0, 2, 0, 3, 0, 1, 0, 5, 
 9, 0, 0, 0, 5, 0, 0, 0, 7, 
 0, 2, 0, 0, 1, 5, 0, 3, 0, 
 0, 0, 0, 6, 0, 0, 2, 0, 0, 
 3, 0, 9, 0, 0, 0, 0, 0, 1), 9, 9)
 puzzle
 sudoku(puzzle)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
