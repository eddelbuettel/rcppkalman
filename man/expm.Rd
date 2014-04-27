\name{expm}
\alias{expm}
\title{Compute the exponential of a matrix}
\usage{
expm(x)
}
\arguments{
  \item{x}{An numeric matrix}
}
\value{
A numeric matrix
}
\description{
This function computes the exponential of a matrix.
}
\details{
This functions calls the \code{expm} function from the
eponymous package \pkg{expm}. This is implemented via a
registered function call, and does not required explicit
linking at the C level. However, the \pkg{expm} package is
imported in order to access its registered function at the
C level.

As the documentation of package \pkg{expm} states, the
underlying implementation borrows from the \pkg{Matrix}
package which itself takes it from GNU Octave.
}
\examples{
## example is from the vignette in package expm
M <- matrix(c(4, 1, 1, 2, 4, 1, 0, 1, 4), 3, 3)

## expected output
expM <- matrix(c(147.8666, 127.7811, 127.7811, 183.7651, 183.7651, 163.6796, 71.79703, 91.88257, 111.96811), 3, 3)

## we only have the expected result to about six digits
all.equal(expm(M), expM, tolerance=1.0e-6)
}
\author{
Dirk Eddelbuettel
}
\seealso{
The \pkg{expm} package and its documentation.
}

