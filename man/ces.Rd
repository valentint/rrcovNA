\name{ces}
\alias{ces}
\docType{data}
\title{ Consumer Expenditure Survey Data }
\description{
    This data set has been derived from the Quarterly Interview Survey of the
    Consumer Expenditure Survey (CES) undertaken by the U.S. Department of Labor,
    Bureau of Labor Statistics and is available at
    \url{https://www.bls.gov/cex/} where also more details about
    this survey can be found. The original data set comprises 869 households in 34
    variables of which one is unique ID, five characterize the size of the
    household, further 6 variables contain other characteristics of the household
    like age, education ethnicity, etc. and 22 variables represent the household
    expenditures. We will consider a reduced set of only 8 expendature variables.
    This reduced data set was analyzed by Hubert at al. (2009)in the context of
    PCA and the first step
    of the analysis showed that all variables are highly skewed.
    They applied the robust PCA method of Serneels and Verdonck based
    on the EM algorithm, since some of the data are incomplete.
}
\usage{data(ces)}
\format{
    A data frame with 869 observations on the following 8 variables:
\describe{
        \item{\code{EXP}}{Total household expenditure}
        \item{\code{FDHO}}{Food and nonalcoholic beverages consumed at home}
        \item{\code{FDAW}}{Food and nonalcoholic beverages consumed away from home}
        \item{\code{SHEL}}{Housing expenditure}
        \item{\code{TELE}}{Telephone services}
        \item{\code{CLOT}}{Clothing}
        \item{\code{HEAL}}{Health care}
        \item{\code{ENT}}{Entertainment}
    }
}
%%\details{}
\source{
     \url{https://www.bls.gov/cex/}
}

\references{
    Hubert, M, Rousseeuw, P.J. and Verdonck, T., (2009). Robust PCA for skewed data and its outlier map,
    \emph{Computational Statistics & Data Analysis}, \bold{53}, 6, pp. 2264-2274
}

\examples{
data(ces)
summary(ces)
plot(ces)
}
\keyword{datasets}
