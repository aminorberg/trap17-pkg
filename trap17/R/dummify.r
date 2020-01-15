#' Make dummy/indicator variables
#'
#' Transform numeric or character variables with multiple levels into dummy/indicators.
#' @param dat Matrix or data frame containing the data
#' @param dumvars Character vector containing the names of the variables to be modified

dummify <- function(dat,
                    dumvars) {

    dums <- c()
    for (i in 1:length(dumvars)) {
        dums <- cbind(dums, 
                      dummies:::dummy(x = dumvars[i], data = dat, drop = FALSE)[, -1])
    }
    return(dums)
}
