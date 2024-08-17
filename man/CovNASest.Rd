\encoding{latin1}
\name{CovNASest}
\alias{CovNASest}
\title{ S Estimates of Multivariate Location and Scatter for incomplete data}
\concept{High breakdown point}
\description{
    Computes S-Estimates of multivariate location and scatter based on Tukey's
    biweight function for incomplete data using a fast algorithm similar to the one proposed by
    Salibian-Barrera and Yohai (2006) for the case of regression. Alternativley, the
    Ruppert's SURREAL algorithm, bisquare or Rocke type estimation can be used.
}
\usage{
    CovNASest(x, bdp = 0.5, arp = 0.1, eps = 1e-5, maxiter = 120,
        nsamp = 500, impMeth = c("norm" , "seq", "rseq"), seed = NULL, 
        trace = FALSE, tolSolve = 1e-13,
        scalefn,
        method = c("sfast", "surreal", "bisquare", "rocke", "suser", "sdet"), control,
        t0, S0, initcontrol)
}
\arguments{
  \item{x}{ a matrix or data frame. }
  \item{bdp}{a numeric value specifying the required
    breakdown point.  Allowed values are between
    \code{(n - p)/(2 * n)} and 1 and the default is \code{bdp=0.5}.
  }
  \item{arp}{a numeric value specifying the asympthotic
    rejection point (for the Rocke type S estimates),
    i.e. the fraction of points receiving zero
    weight (see Rocke (1996)).  Default is \code{arp=0.1}.
  }
  \item{eps}{a numeric value specifying the
    relative precision of the solution of the S-estimate
    (bisquare and Rocke type). Default is to \code{eps=1e-5}.
  }
  \item{maxiter}{maximum number of iterations allowed
    in the computation of the S-estimate (bisquare and Rocke type).
    Default is \code{maxiter=120}.
  }
  \item{nsamp}{the number of random subsets considered. Default is \code{nsamp = 500}.}
  \item{impMeth}{ select imputation method to use - choose one of "norm" , "seq" or "rseq". The default is "norm"}   
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.}
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in
    \code{\link{mahalanobis}}.}
  \item{scalefn}{\code{\link{function}} to compute a robust scale
    estimate or character string specifying a rule determining such a
    function. Used for computing the "deterministic" S-estimates (\code{method="sdet"}).
    If \code{scalefn} is missing or is \code{NULL}, 
    the function is selected depending on the data 
    set size, following the recomendation of Hubert et al. (2012) - 
    \code{\link[robustbase]{Qn}} if \code{n <= 1000} and \code{\link[robustbase]{scaleTau2}} otherwise.}
  \item{method}{ Which algorithm to use: 'sfast'=FAST-S, 'surreal'=SURREAL, 'bisquare', 
  'rocke' or 'sdet', which will invoke the deterministic algorihm of Hubert et al. (2012).}
  \item{control}{ a control object (S4) of class \code{\link[rrcov]{CovControlSest-class}}
    containing estimation options - same as these provided in the fucntion
    specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
  \item{t0}{ optional initial HBDP estimate for the center }
  \item{S0}{ optional initial HBDP estimate for the covariance matrix }
  \item{initcontrol}{ optional control object to be used for computing the initial HBDP estimates }
}
\details{
    Computes biweight multivariate S-estimator of location and scatter. The computation will be
    performed by one of the following algorithms:
    \describe{
        \item{FAST-S}{An algorithm similar to the one proposed by Salibian-Barrera and Yohai (2006) for the case of regression}
        \item{SURREAL}{Ruppert's SURREAL algorithm when \code{method} is set to 'surreal'}
        \item{BISQUARE}{Bisquare S-Estimate with \code{method} set to 'bisquare'}
        \item{ROCKE}{Rocke type S-Estimate with \code{method} set to 'rocke'}.
    }
}
\value{
  An S4 object of class \code{\linkS4class{CovNASest}} which is a subclass of the
  virtual class \code{\linkS4class{CovNARobust}}.
}
\references{
  M. Salibian-Barrera and V. Yohai (2006) A fast algorithm for S-regression
  estimates, \emph{Journal of Computational and Graphical Statistics}, \bold{15},
  414--427.

  R. A. Maronna, D. Martin and V. Yohai (2006). \emph{Robust Statistics: Theory and Methods}.
  Wiley, New York.

  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  <doi:10.18637/jss.v032.i03>.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at},
Matias Salibian-Barrera \email{matias@stat.ubc.ca} and
Victor Yohai \email{vyohai@dm.uba.ar}.  See also the code from
Kristel Joossens, K.U. Leuven, Belgium and Ella Roelant, Ghent University, Belgium.
}
%\note{}
%\seealso{}
\examples{

library(rrcov)
data(bush10)
CovNASest(bush10)

## the following four statements are equivalent
c0 <- CovNASest(bush10)
c1 <- CovNASest(bush10, bdp = 0.25)
c2 <- CovNASest(bush10, control = CovControlSest(bdp = 0.25))
c3 <- CovNASest(bush10, control = new("CovControlSest", bdp = 0.25))

## direct specification overrides control one:
c4 <- CovNASest(bush10, bdp = 0.40,
             control = CovControlSest(bdp = 0.25))
c1
summary(c1)

## Use the SURREAL algorithm of Ruppert
cr <- CovNASest(bush10, method="surreal")
cr

## Use Bisquare estimation
cr <- CovNASest(bush10, method="bisquare")
cr

## Use Rocke type estimation
cr <- CovNASest(bush10, method="rocke")
cr


}
\keyword{robust}
\keyword{multivariate}
