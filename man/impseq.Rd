\name{impSeq}
\alias{impSeq}
\title{
Sequential imputation of missing values
}
\description{
Impute missing multivariate data using sequential algorithm
}
\usage{
impSeq(x)
}
\arguments{
\item{x}{
the original incomplete data matrix.
}}
\value{
a matrix of the same form as \code{x}, but with all missing values filled in sequentially.
}
\details{
\code{SEQimpute} starts from a complete subset of the data set \code{Xc} and estimates 
sequentially the missing values in an incomplete observation, 
say x*, by minimizing the determinant of the covariance of the augmented 
data matrix X* = [Xc; x']. Then the observation x* is added to the complete data matrix
and the algorithm continues with the next observation with missing values.
}
\references{
    S. Verboven, K. Vanden Branden and P. Goos (2007).
    Sequential imputation for missing values.
    \emph{Computational Biology and Chemistry}, bold{31},  320--327.

}
%\seealso{}
\examples{
    data(bush10)
    impSeq(bush10) # impute squentially missing data
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

