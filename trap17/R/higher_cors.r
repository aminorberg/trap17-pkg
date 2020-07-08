#' Higher-level correlations in predictions
#' Compute genoptype-population-level correlations in predictions
#' @param dat Original data
#' @param preds Predicted occurrences
#' @export

higher_cors <- function(dat,
                        preds) {
                        
    form <- formula(cbind(Clo, 
                          Be, 
                          Cap, 
                          Cau, 
                          En) ~ Population + Genotype)

    pooled_occ_true <- as.matrix(aggregate(formula = form, 
                                           data = cbind(dat$X_pooled, dat$Y_pooled), 
                                           FUN = sum))[, -c(1:2)]
    corrs <- matrix(NA, 
                    nrow = dim(preds)[3],
                    ncol = ncol(dat$Y_pooled))
    for (j in 1:dim(preds)[3]) {
        tmp <- preds[, , j]
        colnames(tmp) <- colnames(dat$Y_pooled)
        tmp <- cbind(tmp, dat$X_pooled)
        pred_cors <- as.matrix(aggregate(formula = form, 
                               data = tmp, 
                               FUN = sum))[, -c(1:2)]
        corrs[j,] <- diag(apply(pred_cors, 2, cor, pooled_occ_true, method = "spearman"))
    }
    return(corrs)
}

