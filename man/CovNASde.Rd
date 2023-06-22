\encoding{latin1}
\name{CovNASde}
\Rdversion{1.1}
\alias{CovNASde}
%
\title{Stahel-Donoho Estimates of Multivariate Location and Scatter for incomplete data}
\description{
    Compute a robust estimate of location and scale using the Stahel-Donoho projection based estimator
}
\usage{
CovNASde(x, nsamp, maxres, tune = 0.95, eps = 0.5, prob = 0.99, 
impMeth = c("norm" , "seq", "rseq"), seed = NULL, trace = FALSE, control)
}
\arguments{
  \item{x}{ a matrix or data frame.
}
  \item{nsamp}{a positive integer giving the number of resamples required;
    \code{nsamp} may not be reached if too many of the \code{p}-subsamples,
    chosen out of the observed vectors, are in a hyperplane.
    If \code{nsamp = 0} all possible subsamples are taken.
    If \code{nsamp} is omitted, it is calculated to provide a breakdown point
    of \code{eps} with probability \code{prob}.}
  \item{maxres}{a positive integer specifying the maximum number of
    resamples to be performed including those that are discarded due to linearly
    dependent subsamples.  If \code{maxres} is omitted it will be set to 2 times \code{nsamp}.}
  \item{tune}{a numeric value between 0 and 1 giving the fraction of the data to receive non-zero weight.
    Defaults to \code{0.95}}
  \item{prob}{a numeric value between 0 and 1 specifying the probability of high breakdown point;
    used to compute \code{nsamp} when \code{nsamp} is omitted. Defaults to \code{0.99}.}
  \item{impMeth}{ select imputation method to use - choose one of "norm" , "seq" or "rseq". The default is "norm"}   
  \item{eps}{a numeric value between 0 and 0.5 specifying the breakdown point; used to compute
    \code{nsamp} when \code{nresamp} is omitted. Defaults to \code{0.5}. }
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.  }
  \item{control}{ a control object (S4) of class \code{\link{CovControlSde-class}}
    containing estimation options - same as these provided in the fucntion
    specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
}

%\details{}
\value{
  An S4 object of class \code{\link{CovNASde-class}} which is a subclass of the
  virtual class \code{\link{CovNARobust-class}}.
}

\references{
  R. A. Maronna and V.J. Yohai (1995) The Behavior of the Stahel-Donoho Robust Multivariate
  Estimator. \emph{Journal of the American Statistical Association} \bold{90} (429), 330--341.

  R. A. Maronna, D. Martin and V. Yohai (2006). \emph{Robust Statistics: Theory and Methods}.
  Wiley, New York.

  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  <doi:10.18637/jss.v032.i03>.
}

\author{
Valentin Todorov \email{valentin.todorov@chello.at}
}

\examples{
data(bush10)
CovNASde(bush10)

## the following four statements are equivalent
c0 <- CovNASde(bush10)
c1 <- CovNASde(bush10, nsamp=2000)
c2 <- CovNASde(bush10, control = CovControlSde(nsamp=2000))
c3 <- CovNASde(bush10, control = new("CovControlSde", nsamp=2000))

## direct specification overrides control one:
c4 <- CovNASde(bush10, nsamp=100,
             control = CovControlSde(nsamp=2000))
c1
summary(c1)
}
\keyword{robust}
\keyword{multivariate}
