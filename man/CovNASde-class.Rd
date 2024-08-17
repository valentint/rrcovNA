\name{CovNASde-class}
\Rdversion{1.1}
\docType{class}
\alias{CovNASde-class}

\title{
Stahel-Donoho Estimates of Multivariate Location and Scatter for incomplete data
}
\description{ This class, derived from the virtual class \code{"CovRobust"}
    accomodates Stahel-Donoho estimates of multivariate location and scatter.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovSde", ...)},
but the usual way of creating \code{CovSde} objects is a call to the function
\code{CovSde} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{iter}, \code{crit}, \code{wt}:}{from the
            \code{"\link[rrcov]{CovRobust-class}"} class.
        }
        \item{\code{call}, \code{cov}, \code{center},
            \code{n.obs}, \code{mah}, \code{method},
            \code{singularity}, \code{X}:}{from the \code{"\link[rrcov]{Cov-class}"} class.
        }
    }
}
\section{Extends}{
    Class \code{"\link[rrcov]{CovRobust-class}"}, directly.
    Class \code{"\link[rrcov]{Cov-class}"}, by class \code{"\link[rrcov]{CovRobust-class}"}.
}
\section{Methods}{
No methods defined with class "CovNASde" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  <doi:10.18637/jss.v032.i03>.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
   \code{\link[rrcov]{CovSde}}, \code{\link[rrcov]{Cov-class}}, \code{\link[rrcov]{CovRobust-class}}
}
\examples{
showClass("CovNASde")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
