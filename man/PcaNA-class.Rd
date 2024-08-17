\name{PcaNA-class}
\Rdversion{1.1}
\docType{class}
\alias{PcaNA-class}
\alias{getQuan,PcaNA-method}

\title{Class "PcaNA" Principal Components for incomplete data }
\description{
Contains the results of the computations of  classical and robust principal components for incomplete data using an EM algorithm as descibed by Serneels and Verdonck (2008)}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PcaNA", ...)} but the
    usual way of creating \code{PcaNA} objects is a call to the function
    \code{PcaNA} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{call}, \code{center}, \code{scale}, \code{loadings},
            \code{eigenvalues}, \code{scores}, \code{k},
            \code{sd}, \code{od}, \code{cutoff.sd}, \code{cutoff.od},
            \code{flag}, \code{n.obs}:}{
            from the \code{"\link[rrcov]{Pca-class}"} class.
        }
        \item{\code{Ximp}:}{the data matrix with imputed missing values}
    }
}

\section{Extends}{
Class \code{"\link[rrcov]{Pca-class}"}, directly.
}
\section{Methods}{
  \describe{
    \item{getQuan}{\code{signature(obj = "PcaNA")}: ... }
     }
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
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link[rrcov]{PcaRobust-class}}, \code{\link[rrcov]{Pca-class}}, 
     \code{\link[rrcov]{PcaClassic}}, \code{\link[rrcov]{PcaClassic-class}}
}

\examples{
showClass("PcaNA")
}
\keyword{robust}
\keyword{multivariate}
