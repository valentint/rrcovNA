\name{impSeq}
\alias{impSeq}
\title{
Sequential imputation of missing values
}
\description{
Impute missing multivariate data using sequential algorithm
}
\usage{
impSeq(x, norm_impute=FALSE, check_data=FALSE, verbose=TRUE)
}
\arguments{
\item{x}{the original incomplete data matrix.}
\item{norm_impute}{If there are not enough complete observations and \code{norm_impute=TRUE}
    the data matrix will be imputed using the multivariate normal model,
    the outlyingness will be computed for the completed matrix and the result will be returned.
    If \code{norm_impute=FALSE} (the default) the algorithm will try to impute 
    sufficient number of observations using the multivariate normal model and then will
    start the sequential imputation.}
\item{check_data}{whether to check the variables: only numeric, non discrete, 
    with less than 50\% NAs and with non-zero MAD. The default is \code{check_data=FALSE},
    meaning that no checks will be performed.}
\item{verbose}{whether to write messages about the checking of the data. By default 
    \code{verbose=TRUE}, i.e. messages will be written. }

}
\value{
A list containing the following elements:
    \item{\code{x}}{a matrix of the same form as \code{x}, but with all missing values filled in sequentially.}
    \item{\code{colInAnalysis}}{the column indices of the columns used in the analysis.}
    \item{\code{namesNotNumeric}}{the names of the variables which are not numeric.}
    \item{\code{namesNAcol}}{names of the columns left out due to too many NA's.}
    \item{\code{namesDiscrete}}{names of the discrete variables.}
    \item{\code{namesZeroScale}}{names of the variables with zero scale.}
}
\details{
\code{SEQimpute} starts from a complete subset of the data set \code{Xc} and estimates 
sequentially the missing values in an incomplete observation, 
say x*, by minimizing the determinant of the covariance of the augmented 
data matrix X* = [Xc; x']. Then the observation x* is added to the complete data matrix
and the algorithm continues with the next observation with missing values.
}
\references{
    S. Verboven, K. Vanden Branden and P. Goos (2007).
    Sequential imputation for missing values.
    \emph{Computational Biology and Chemistry}, \bold{31},  320--327.

}
%\seealso{}
\examples{
    data(bush10)
    impSeq(bush10) # impute squentially missing data
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

