\name{impNorm}
\alias{impNorm}
\title{
Impute missing multivariate normal data
}
\description{
Draws missing elements of a data matrix under the multivariate normal
model and a user-supplied parameter
}
\usage{
impNorm(x)
}
\arguments{
\item{x}{
the original incomplete data matrix.
}}
\value{
a matrix of the same form as \code{x}, but with all missing values filled in
with simulated values drawn from their predictive distribution given
the observed data and the specified parameter.
}
\details{
This function simply uses \code{imp.norm} from package \code{norm}.
}
\references{
See Section 5.4.1 of Schafer (1996).
}
\seealso{
\code{\link[norm]{prelim.norm}}, \code{\link[norm]{makeparam.norm}}, and \code{\link[norm]{rngseed}}.
}
\examples{
data(bush10)
impNorm(bush10) #impute missing data under the MLE

}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

