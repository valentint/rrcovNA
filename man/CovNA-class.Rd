\name{CovNA-class}
\docType{class}
\alias{CovNA-class}
\alias{getDistance,CovNA-method}
\alias{getFlag,CovNA-method}
\alias{summary,CovNA-method}
\title{Class "CovNA" -- a base class for estimates of multivariate
    location and scatter for incomplete data}
\description{The class \code{CovNA} represents an estimate of the
    multivariate location and scatter of a data set. The objects of class \code{CovNA}
    contain the classical estimates and serve as base for deriving other
    estimates, i.e. different types of robust estimates.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovNA", ...)},
    but the usual way of creating \code{CovNA} objects is a call to the function
    \code{CovNA} which serves as a constructor.
}
\section{Slots}{
    \describe{
    \item{\code{call}:}{Object of class \code{"language"} }
    \item{\code{cov}:}{covariance matrix}
    \item{\code{center}:}{location}
    \item{\code{n.obs}:}{number of observations used for the computation of the estimates}
    \item{\code{mah}:}{mahalanobis distances}
    \item{\code{det}:}{determinant}
    \item{\code{flag}:}{flags (FALSE if suspected an outlier)}
    \item{\code{method}:}{a character string describing the method used to
        compute the estimate: "Classic"}
    \item{\code{singularity}:}{a list with singularity information for the
        covariance matrix (or \code{NULL} of not singular)}
    \item{\code{X}:}{data}
  }
}
\section{Extends}{
    Class \code{"\linkS4class{Cov}"}, directly.
}
\section{Methods}{
  \describe{
    \item{getDistance}{\code{signature(obj = "CovNA")}: distances }
    \item{getFlag}{\code{signature(obj = "CovNA")}: Flags observations as outliers if the corresponding mahalanobis distance is larger then \code{qchisq(prob, p)} where \code{prob} defaults to 0.975. }
    \item{summary}{\code{signature(object = "CovNA")}: calculate summary information }
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
showClass("CovNA")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
