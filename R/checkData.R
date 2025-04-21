checkData <- function(X, check_data=TRUE, fracNA=0.5, numDiscrete=3, precScale=1e-12, verbose=TRUE) {
## This function checks the data X, and removes
## columns and rows that do not satisfy the conditions.
##
## fracNA      : only consider columns and rows with fewer NAs than this.
## numDiscrete : a column that takes on numDiscrete or fewer values
##               will be considered discrete and not used in the analysis.
## precScale   : only consider columns whose scale is > precScale.
##               Here scale is measured by the median absolute deviation.

  
    ## Add column names and row names if these were not given:
    if(is.matrix(X)) 
        X <- data.frame(X)
    ## This keeps the row names and column names if they exist, else creates
    ## them as 1, 2, 3, ... for rows and X1, X2, X3,... for columns.

    x <- X                                  # keep the original variables in X
    rowInAnalysis <- 1:nrow(X)              # all rows
    colInAnalysis <- which(sapply(x, is.numeric))  # only numeric columns
    namesNotNumeric <- NULL  
    namesNAcol <- NULL
    namesDiscrete <- NULL
    namesZeroScale <- NULL
    
    numgoodcol <- length(colInAnalysis)
    vecNotNumeric <- which(!(colInAnalysis %in% 1:ncol(X)))
    numbadcol <- length(vecNotNumeric)      # can be 0

    if(numbadcol > 0) {
        namesNotNumeric <- colnames(x)[vecNotNumeric]
        x <- X[, colInAnalysis]
        if(verbose) {
            cat("\nThe data contained ", numbadcol, " non-numeric columns:\n")
            print(namesNotNumeric)   
        }
        if(numgoodcol == 0)
            stop("No columns remain, cannot continue!")
    } 
    
    if(check_data) {
        ##  2. Exclude variables with more than 50% NAs
        indNAcol <- which(apply(x, 2, function(x) length(which(is.na(x))))/nrow(x) >= fracNA)
        if(length(indNAcol) > 0) {
            namesNAcol <- colnames(x)[indNAcol]
            x <- x[, -indNAcol]
            ngoodcol <- ncol(x)
            if(verbose) {
                cat("\nThe data contained",  length(indNAcol), "column(s) with over", 100*fracNA, "% of NAs:\n")
                print(namesNAcol)
            }
            if(ngoodcol <= 0)
                stop("No columns remain, cannot continue!")
        }
    
        ## 3. Exclude discrete columns with 3 or fewer values
        indDiscrete <- which(apply(x, 2, function(x) sum(!is.na(unique(x)))) <= numDiscrete)
        if(length(indDiscrete) > 0) {
            namesDiscrete <- colnames(x)[indDiscrete]
            x <- x[, -indDiscrete]
            ngoodcol <- ncol(x)
            if(verbose) {
                cat("\nThe data contained",  length(indDiscrete), "discrete columns with", numDiscrete, "or fewer values:\n")
                print(namesDiscrete)
            }
            if(ngoodcol <= 0)
                stop("No columns remain, cannot continue!")
        }
    
        ## 4. Exclude columns for which the median absolute deviation is
        ##    zero. This is equivalent to saying that 50% or more of its
        ##    values are equal.    
        indZeroScale <- which(apply(x, 2, mad, na.rm=TRUE) <= precScale)
        if(length(indZeroScale) > 0) {
            namesZeroScale <- colnames(x)[indZeroScale]
            x <- x[, -indZeroScale]
            ngoodcol <- ncol(x)
            if(verbose) {
                cat("\nThe data contained",  length(indZeroScale), "columns with zero or tiny median absolute deviation:\n")
                print(namesZeroScale)
            }
            if(ngoodcol <= 0)
                stop("No columns remain, cannot continue!")
        }
    }    
    
    colInAnalysis <- which(colnames(X) %in% colnames(x))
    x <- X[rowInAnalysis, colInAnalysis]

    return(list(colInAnalysis=colInAnalysis,
##              rowInAnalysis=rowInAnalysis,
              namesNotNumeric=namesNotNumeric,
##              namesCaseNumber = namesCaseNumber,
              namesNAcol=namesNAcol,
##              namesNArow = namesNArow,
              namesDiscrete=namesDiscrete,
              namesZeroScale=namesZeroScale))
}    
