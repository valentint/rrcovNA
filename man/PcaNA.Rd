\name{PcaNA}
\alias{PcaNA}
\alias{PcaNA.formula}
\alias{PcaNA.default}

\title{ Classical or robust Principal Components for incomplete data}
\description{
Computes classical and robust principal components for incomplete data using an EM algorithm as descibed by Serneels and Verdonck (2008)
}
\usage{
PcaNA(x, ...)
\method{PcaNA}{default}(x, k = ncol(x), kmax = ncol(x), conv=1e-10, maxiter=100, 
    method=c("cov", "locantore", "hubert", "grid", "proj", "class"), cov.control=NULL,
    scale = FALSE, signflip = TRUE, crit.pca.distances = 0.975, trace=FALSE, \dots)
\method{PcaNA}{formula}(formula, data = NULL, subset, na.action, \dots)
}
\arguments{
  \item{formula}{a formula with no response variable, referring only to
    numeric variables.}
  \item{data}{an optional data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.}
  \item{subset}{an optional vector used to select rows (observations) of the
    data matrix \code{x}.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset. The default is \code{\link{na.omit}}.}
  \item{\dots}{arguments passed to or from other methods.}
  \item{x}{a numeric matrix (or data frame) which provides
    the data for the principal components analysis.}
  \item{k}{number of principal components to compute. If \code{k} is missing, 
    or \code{k = 0}, the algorithm itself will determine the number of 
    components by finding such \code{k} that \eqn{l_k/l_1 >= 10.E-3} and 
    \eqn{\Sigma_{j=1}^k l_j/\Sigma_{j=1}^r l_j >= 0.8}. 
    It is preferable to investigate the scree plot in order to choose the number 
    of components and then run again. Default is \code{k=ncol(x)}. }
  \item{kmax}{maximal number of principal components to compute.
    Default is \code{kmax=10}. If \code{k} is provided, \code{kmax} 
    does not need to be specified, unless \code{k} is larger than 10.}
  \item{conv}{convergence criterion for the EM algorithm.
    Default is \code{conv=1e-10}.}
  \item{maxiter}{maximal number of iterations for the EM algorithm.
    Default is \code{maxiter=100}.}
  \item{method}{which PC method to use (classical or robust) - "class" means classical PCA 
    and one of the following "locantore", "hubert", "grid", "proj", "cov" specifies a 
    robust PCA method. If the method is "cov" - i.e. PCA based on a robust covariance matrix - 
    the argument \code{cov.control} can specify which method for computing the 
    (robust) covariance matrix will be used. 
    Default is \code{method="locantore"}.}
  \item{cov.control}{control object in case of robust PCA based on a robust covariance matrix.}
  \item{scale}{a logical value indicating whether the variables should be 
  scaled to have unit variance (only possible if there are no constant 
  variables). As a scale function \code{mad} is used but alternatively, a vector of length equal 
  the number of columns of x can be supplied. The value is passed to 
  scale and the result of the scaling is stored in the \code{scale} slot. 
  Default is \code{scale = FALSE}}
  \item{signflip}{a logical value indicating wheather to try to solve the sign indeterminancy of the loadings -   
  ad hoc approach setting the maximum element in a singular vector to be positive. Default is \code{signflip = FALSE}}
  \item{crit.pca.distances}{criterion to use for computing the cutoff values for the orthogonal and score distances. Default is 0.975.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}
\details{
\code{PcaNA}, serving as a constructor for objects of class \code{\linkS4class{PcaNA}} 
is a generic function with "formula" and "default" methods. For details see the relevant references.
}
\value{
  An S4 object of class \code{\linkS4class{PcaNA}} which is a subclass of the 
  virtual class \code{\link[rrcov]{Pca-class}}. 
}
\references{

  Serneels S & Verdonck T (2008),
  Principal component analysis for data containing outliers and missing elements.
  \emph{Computational Statistics and Data Analisys}, \bold{52}(3), 1712--1727 .

  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  <doi:10.18637/jss.v032.i03>.

}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
%\note{}
\examples{
## 1. With complete data
## PCA of the bushfire data
    data(bushfire)
    pca <- PcaNA(bushfire)
    pca

## Compare with the classical PCA
    prcomp(bushfire)

## or  
    PcaNA(bushfire, method="class")
    
## If you want to print the scores too, use
    print(pca, print.x=TRUE)

## Using the formula interface
    PcaNA(~., data=bushfire)

## To plot the results:

    plot(pca)                    # distance plot
    pca2 <- PcaNA(bushfire, k=2)  
    plot(pca2)                   # PCA diagnostic plot (or outlier map)
    
## Use the standard plots available for for prcomp and princomp
    screeplot(pca)    
    biplot(pca)  

################################################################      
## 2. Now the same wit incomplete data - bush10
    data(bush10)
    pca <- PcaNA(bush10)
    pca

## Compare with the classical PCA
    PcaNA(bush10, method="class")
    
## If you want to print the scores too, use
    print(pca, print.x=TRUE)

## Using the formula interface
    PcaNA(~., data=as.data.frame(bush10))

## To plot the results:

    plot(pca)                    # distance plot
    pca2 <- PcaNA(bush10, k=2)  
    plot(pca2)                   # PCA diagnostic plot (or outlier map)
    
## Use the standard plots available for for prcomp and princomp
    screeplot(pca)    
    biplot(pca)    
    
}
\keyword{robust}
\keyword{multivariate}
