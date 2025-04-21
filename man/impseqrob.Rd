\name{impSeqRob}
\alias{impSeqRob}
\title{
Robust sequential imputation of missing values
}
\description{
Impute missing multivariate data using robust sequential algorithm
}
\usage{
impSeqRob(x, alpha=0.9, norm_impute=FALSE, check_data=FALSE, verbose=TRUE)
}
\arguments{
\item{x}{the original incomplete data matrix.}
\item{alpha}{the number of regular genes, \code{(1-alpha)} measures the 
    fraction of outliers the algorithm should resist. Any value between 0.5 
    and 1 may be specified. The default is \code{alpha=0.9}.}
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
    \item{\code{outl}}{outlyingness computed for all observations.}
    \item{\code{flag}}{flag of outliers.}
    \item{\code{colInAnalysis}}{the column indices of the columns used in the analysis.}
    \item{\code{namesNotNumeric}}{the names of the variables which are not numeric.}
    \item{\code{namesNAcol}}{names of the columns left out due to too many NA's.}
    \item{\code{namesDiscrete}}{names of the discrete variables.}
    \item{\code{namesZeroScale}}{names of the variables with zero scale.}
}
\details{
The nonrobust version \code{SEQimpute} starts from a complete subset of the data set \code{Xc} and estimates 
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
    \emph{Computational Biology and Chemistry}, \bold{31},  320--327.
    
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

