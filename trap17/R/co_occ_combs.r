#' CALCULATE OCCURRENCE COMBINATIONS
#'
#' Calculate combinations of species occurrences
#' @param partition A matrix giving the variables from which to form the combinations
#' @param Y A list of matrices of (species) occurrences
#' @export

co_occ_combs <- function(partition,
                         Y_arr) 
{

    Xpartition <- as.data.frame(partition)
    for (i in 1:ncol(Xpartition)) {
        Xpartition[,i] <- paste0(colnames(Xpartition)[i], Xpartition[,i])    
    }

    calc_combs <- function (Y, Xpart) {        
        combpaste <- function(y, sp_names) {
             tmp1 <- paste(sp_names[as.logical(y)], 
                           collapse = "+")
            return(tmp1)           
        }
        y_df <- as.data.frame(Y)           
        tmp11 <- apply(y_df, MARGIN = 1, FUN = combpaste, sp_names = dimnames(y_df)[[2]])
        Combinations <- matrix(tmp11, ncol = 1)

        res_iter <- switch(ncol(Xpart), 
                           "1" = cbind(Combinations, Xpart[,1]), 
                           "2" = cbind(Combinations, Xpart[,1], Xpart[,2]), 
                           "3" = cbind(Combinations, Xpart[,1], Xpart[,2], Xpart[,3]),
                           "4" = cbind(Combinations, Xpart[,1], Xpart[,2], Xpart[,4]),
                           "5" = cbind(Combinations, Xpart[,1], Xpart[,2], Xpart[,3], Xpart[,4], Xpart[,5]))
        return(res_iter)
    }

    tmp_res <- apply(Y_arr, MARGIN = 3, FUN = calc_combs, Xpart = Xpartition)
    
    
    res_iters <- NA
    for (i in 1:ncol(tmp_res)) {
        tmp1 <- matrix(tmp_res[,i], ncol = (ncol(Xpartition)+1), nrow = (length(tmp_res[,i])/3))
        res_iters <- rbind(res_iters, tmp1)
    }

    res_iters <- res_iters[-1,]
    res_iters[which(res_iters == "")] <- "Empty"

    res <- switch(ncol(Xpartition), 
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
