\name{SummaryCovNARobust-class}
\docType{class}
\alias{SummaryCovNARobust-class}
\alias{show,SummaryCovNARobust-method}

\title{Class "SummaryCovNARobust" - summary of "CovNARobust" objects }
\description{    Summary information for CovRobust objects meants for printing by 'show' }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SummaryCovNARobust", ...)}, 
but most often by invoking 'summary' on an "CovNA" object. They contain values 
meant for printing by 'show'.
}
\section{Slots}{
No Slots defined with class \code{"SummaryCovNARobust"} in the signature.
}
\section{Extends}{
Class \code{"SummaryCovNA"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SummaryCovNARobust")}: ... }
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
  \code{\link[rrcov]{CovRobust-class}},   \code{\link[rrcov]{SummaryCov-class}}
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
cv <- CovMest(hbk.x)
cv
summary(cv)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
