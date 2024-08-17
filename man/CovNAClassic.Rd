\name{CovNAClassic}
\alias{CovNAClassic}
\alias{CovNA}
\title{Classical Estimates of Multivariate Location and Scatter for incomplete data (EM Algorithm)}
\concept{High breakdown point}
\description{
    Computes the classical estimates of multivariate location and scatter.
    Returns an S4 class \code{CovNAClassic} with the estimated \code{center},
    \code{cov}, Mahalanobis distances and weights based on these distances.
}
\usage{
    CovNAClassic(x, unbiased=TRUE)
    CovNA(x, unbiased=TRUE)
}
\arguments{
  \item{x}{a matrix or data frame. As usual, rows are observations and
        columns are variables. }
  \item{unbiased}{ whether to return the unbiased estimate of
    the covariance matrix. Default is  \code{unbiased = TRUE}}
}
%\details{}
\value{
  An object of class \code{"CovNAClassic"}.
}
\references{
  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  <doi:10.18637/jss.v032.i03>.
}
%\note{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\seealso{
  \code{\link[rrcov]{Cov-class}}, \code{\link[rrcov]{CovClassic-class}}, \code{\link{CovNAClassic-class}}
}
\examples{
data(bush10)
cv <- CovNAClassic(bush10)
cv
summary(cv)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
