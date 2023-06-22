\name{bush10}
\alias{bush10}
\docType{data}
\title{ Campbell Bushfire Data with added missing data items (10 percent)}
\description{
  This data set is based on the bushfire data set which was used by
  Campbell (1984) to locate bushfire scars - see  \code{\link[robustbase]{bushfire}}
  in package \code{robustbase}. The original dataset contains satelite
  measurements on five frequency bands, corresponding to each of 38 pixels.
  The data set is very well studied (Maronna and Yohai, 1995; Maronna
and Zamar, 2002). There are 12 clear outliers: 33-38, 32, 7-11 and 12 and 13 are
suspect.
}
\usage{data(bush10)}
\format{
  A data frame with 38 observations on 6 variables.

  The original data set consists of 38 observations in 5 variables.
  Based on it four new data sets are created in which some of the data
  items are replaced by missing values with a simple "missing completely
  at random " mechanism. For this purpose independent Bernoulli trials
  are realized for each data item with a probability of success 0.1
  where success means that the corresponding item is set to missing.)
  %\describe{}
}
%\Note{}
\source{
  Maronna, R.A. and Yohai, V.J. (1995)
  The Behavoiur of the Stahel-Donoho Robust Multivariate Estimator.
  \emph{Journal of the American Statistical Association} \bold{90}, 330--341.

  Beguin, C. and Hulliger, B. (2004) Multivariate outlier detection in incomplete
    survey data: the epidemic algorithm and transformed rank correlations.
    \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)} \bold{127}, 2, 275--294.
}
%\seealso{}
\examples{
## The following code will result in exactly the same output
##  as the one obtained from the original data set
data(bush10)
plot(bush10)
CovNAMcd(bush10)


\dontrun{
##  This is the code with which the missing data were created:
##  Creates a data set with missing values (for testing purposes)
##  from a complete data set 'x'. The probability of
##  each item being missing is 'pr'.
##
getmiss <- function(x, pr=0.1){
    library(Rlab)
    n <- nrow(x)
    p <- ncol(x)
    bt <- rbern(n*p, pr)
    btmat <- matrix(bt, nrow=n)
    btmiss <- ifelse(btmat==1, NA, 0)
    x+btmiss
}
}

}
\keyword{datasets}
