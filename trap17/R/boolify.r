#' Make boolean variables
#'
#' Make numeric or character variables boolean
#' @param dat Matrix or data frame containing the data

boolify <- function(dat,
                    boolvars) {

    bools <- c()
    for (i in 1:length(boolvars)) {
        bools <- cbind(bools, 1 * (dat[, boolvars[i]] > 0))
    }
    bools <- matrix(bools, ncol = length(boolvars))
    colnames(bools) <- boolvars
    return(bools)
}
