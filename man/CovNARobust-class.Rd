\name{CovNARobust-class}
\Rdversion{1.1}
\docType{class}
\alias{CovNARobust-class}
\alias{summary,CovNARobust-method}
\alias{plot,CovNARobust,missing-method}
\title{Class "CovNARobust" - virtual base class for robust
 estimates of multivariate location and scatter for incomplete data}

\description{  \code{CovNARobust} is a virtual base class used for deriving the concrete classes
representing different robust estimates of multivariate location and scatter for incomplete data. Here are implemeted the
standard methods common for all robust estimates like \code{show}, \code{summary} and \code{plot}.
The derived classes can override these methods and can define new ones.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
    \describe{
        \item{\code{iter}:}{number of iterations used to compute the estimates}
        \item{\code{crit}:}{value of the criterion function}
        \item{\code{wt}:}{weights}
        \item{\code{call}, \code{cov}, \code{center},
            \code{n.obs}, \code{mah}, \code{method},
            \code{singularity}, \code{X}:}{from the \code{"\linkS4class{Cov}"} class.
        }
    }
}
\section{Extends}{
    Classes \code{"\linkS4class{CovNA}"} and \code{"\linkS4class{CovRobust}"}, directly.
}
\section{Methods}{
\describe{
    \item{plot}{\code{signature(x = "CovNARobust")}: plot the object }
    \item{summary}{\code{signature(object = "CovNARobust")}: display additional information for the object }
}
}
\references{
  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  <doi:10.18637/jss.v032.i03>.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\note{}

\seealso{
  \code{\link{CovNA-class}}, \code{\link{CovNAMcd-class}}, \code{\link{CovNAOgk-class}}, \code{\link{CovNASde-class}}, , \code{\link{CovNASest-class}}
}
\examples{
     data(hbk)
     hbk.x <- data.matrix(hbk[, 1:3])
     cv <- CovMest(hbk.x)               # it is not possible to create an object of
                                        # class CovRobust, since it is a VIRTUAL class
     cv
     summary(cv)                        # summary method for class CovRobust
     plot(cv)                           # plot method for class CovRobust
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
