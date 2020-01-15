#' CALCULATE OCCURRENCE COMBINATIONS
#'
#' Calculate combinations of species occurrences
#' @param partition A matrix giving the variables from which to form the combinations
#' @param Y A list of matrices of (species) occurrences
#' @export

co_occ_combs <- function(partition,
                         Y_arr) 
{

    X <- as.data.frame(partition)
    for (i in 1:ncol(X)) {
        X[,i] <- paste0(colnames(X)[i], X[,i])    
    }

    res_iters <- NA
    for (j in 1:dim(Y_arr)[3]) {
        Y <- as.data.frame(Y_arr[,,j])   

        Combinations <- paste(dimnames(Y)[2][[1]][as.logical(Y[1,])], 
                        collapse = "+")
        for (k in 2:dim(Y)[1]) {
             tmp1 <- paste(dimnames(Y)[2][[1]][as.logical(Y[k,])], 
                           collapse = "+")
             Combinations <- rbind(Combinations, tmp1)
        }

#        res_iter <- switch(ncol(X), 
#                         "1" = table(Combinations, X[,1]), 
#                         "2" = table(Combinations, X[,1], X[,2]), 
#                         "3" = table(Combinations, X[,1], X[,2], X[,3]),
#                         "4" = table(Combinations, X[,1], X[,2], X[,4]),
#                         "5" = table(Combinations, X[,1], X[,2], X[,3], X[,4], X[,5]))
        tmp2 <- switch(ncol(X), 
                         "1" = cbind(Combinations, X[,1]), 
                         "2" = cbind(Combinations, X[,1], X[,2]), 
                         "3" = cbind(Combinations, X[,1], X[,2], X[,3]),
                         "4" = cbind(Combinations, X[,1], X[,2], X[,4]),
                         "5" = cbind(Combinations, X[,1], X[,2], X[,3], X[,4], X[,5]))
        res_iters <- rbind(res_iters, tmp2)
    }
    res_iters <- res_iters[-1,]
    res_iters[which(res_iters == "")] <- "Empty"

    res <- switch(ncol(X), 
                     "1" = table(res_iters[,1], 
                                 res_iters[,2]), 
                     "2" = table(res_iters[,1], 
                                 res_iters[,2], 
                                 res_iters[,3]), 
                     "3" = table(res_iters[,1], 
                                 res_iters[,2], 
                                 res_iters[,3], 
                                 res_iters[,4]),
                     "4" = table(res_iters[,1], 
                                 res_iters[,2], 
                                 res_iters[,3], 
                                 res_iters[,4], 
                                 res_iters[,5]),
                     "5" = table(res_iters[,1], 
                                 res_iters[,2], 
                                 res_iters[,3], 
                                 res_iters[,4], 
                                 res_iters[,5], 
                                 res_iters[,6]))
                                 
    return(res)
}
