\name{CovNAClassic-class}
\docType{class}
\alias{CovNAClassic-class}
\alias{plot,CovNAClassic,missing-method}
\title{Class "CovNAClassic" - classical estimates of multivariate
    location and scatter for incomplete data (EM algorithm)}
\description{ The class \code{CovNAClassic} represents an estimate of the
    multivariate location and scatter of an incomplete data set. The class \code{CovNAClassic}
    objects contain the classical estimates.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovNAClassic", ...)},
    but the usual way of creating \code{CovNAClassic} objects is a call to the function
    \code{CovNAClassic} which serves as a constructor.
}
\section{Slots}{
    \describe{
    \item{\code{call}:}{Object of class \code{"language"} }
    \item{\code{cov}:}{covariance matrix}
    \item{\code{center}:}{location}
    \item{\code{n.obs}:}{number of observations used for the computation of the estimates}
    \item{\code{mah}:}{mahalanobis distances}
    \item{\code{method}:}{a character string describing the method used to
        compute the estimate: "Classic"}
    \item{\code{singularity}:}{a list with singularity information for the
        ocvariance matrix (or \code{NULL} of not singular)}
    \item{\code{X}:}{data}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "CovNAClassic")}: plot the object }
   }
}
\references{
  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  <doi:10.18637/jss.v032.i03>.
}
%\note{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\seealso{}
\examples{
data(bush10)
cv <- CovNAClassic(bush10)
cv
summary(cv)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
