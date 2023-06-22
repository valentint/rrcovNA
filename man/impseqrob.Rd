\name{impSeqRob}
\alias{impSeqRob}
\title{
Robust sequential imputation of missing values
}
\description{
Impute missing multivariate data using robust sequential algorithm
}
\usage{
impSeqRob(x, alpha=0.9)
}
\arguments{
\item{x}{the original incomplete data matrix.}
\item{alpha}{.The default is \code{alpha=0.9}.}

}
\value{
a matrix of the same form as \code{x}, but with all missing values filled in sequentially.
}
\details{
\code{SEQimpute} starts from a complete subset of the data set \code{Xc} and estimates 
sequentially the missing values in an incomplete observation, 
say x*, by minimizing the determinant of the covariance of the augmented 
data matrix X* = [Xc; x']. Then the observation x* is added to the complete data matrix
and the algorithm continues with the next observation with missing values.
Since \code{SEQimpute} uses the sample mean and covariance matrix it will be vulnerable 
to the influence of outliers and it is improved by plugging in robust estimators of 
location and scatter. One possible solution is to use the outlyingness measure as proposed 
by Stahel (1981) and Donoho (1982) and successfully used for outlier 
identification in Hubert et al. (2005). We can compute the outlyingness measure for 
the complete observations only but once an incomplete observation is imputed (sequentially)  
we could compute the outlyingness measure for it too and use it to decide if this observation 
is an outlier or not. If the outlyingness measure does not exceed a predefined threshold 
the observation is included in the further steps of the algorithm. 
}
\references{
    S. Verboven, K. Vanden Branden and P. Goos (2007).
    Sequential imputation for missing values.
    \emph{Computational Biology and Chemistry}, bold{31},  320--327.
    K. Vanden Branden and S. Verboven (2009).
    Robust Data Imputation.
    \emph{Computational Biology and Chemistry}, \bold{33}, 7--13.
}
%\seealso{}
\examples{
    data(bush10)
    impSeqRob(bush10) # impute squentially missing data
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

