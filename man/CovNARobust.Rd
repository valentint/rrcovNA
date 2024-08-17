\name{CovNARobust}
\alias{CovNARobust}
\title{ Robust Location and Scatter Estimation for incomplete data}
\description{
    Computes a robust multivariate location and scatter estimate with a high
     breakdown point for incomplete data, using one of the available estimators.
}
\usage{
CovNARobust(x, control, impMeth=c("norm" , "seq", "rseq"))
}
\arguments{
  \item{x}{a matrix or data frame. }
  \item{control}{ a control object (S4) for one of the available control classes,
    e.g. \code{\link[rrcov]{CovControlMcd-class}}, \code{\link[rrcov]{CovControlOgk-class}},
    \code{\link[rrcov]{CovControlSest-class}}, etc.,
    containing estimation options. The class of this object defines
    which estimator will be used. Alternatively a character string can be specified
    which names the estimator - one of auto, sde, mcd, ogk, m, mve, sfast, surreal,
    bisquare, rocke. If 'auto' is specified or the argument is missing, the
    function will select the estimator (see below for details)}
  \item{impMeth}{ select imputation method to use - choose one of "norm" , "seq" or "rseq". The default is "norm"}
}
\details{
  This function is based on imputation and than estimation with a selected high breakdown point method. 
  Thus first imputation with the selected method will be performed and then the function \code{CovRobust}will be called. 
  For details see \code{\link[rrcov]{CovRobust}}.
}
\value{
    An object derived from a \code{CovRobust} object, depending on the selected estimator.
}
\references{
  V. Todorov, M. Templ and P. Filzmoser. Detection of multivariate outliers in business
  survey data with incomplete information. \emph{Advances in Data Analysis and Classification}, 
  \bold{5} 37--56, 2011.

  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  <doi:10.18637/jss.v032.i03>.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}
}

\examples{
data(bush10)
CovNARobust(bush10)
CovNARobust(bush10, CovControlSest())
}
\keyword{robust}
\keyword{multivariate}
