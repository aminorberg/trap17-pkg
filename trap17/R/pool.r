#' Pool trap17 sata
#'
#' Pool the trap17 data over the two sampling time points
#' @param dirs Directories
#' @param save_data Save the pooled data
#' @param return_data Return the pooled data

pool <- function(dat, dirs, save_data = FALSE) 

{

    res <- structure(list(Y_pooled = NULL, X_pooled = NULL), 
                     class = "trapdata")

    dfr <- data.frame(cbind(dat$X, dat$Y))

    pool_sum <- aggregate(formula = dat$pool_sum_formula, data = dfr, FUN = sum)
    pool_min <- aggregate(formula = dat$pool_min_formula, data = dfr, FUN = min)
    pool_max <- aggregate(formula = dat$pool_max_formula, data = dfr, FUN = max)
    pooled1 <- merge(pool_sum, pool_min, by = dat$poolvar, all = TRUE)
    pooled <- merge(pooled1, pool_max, by = dat$poolvar, all = TRUE)

    Y_pooled <- as.matrix(pooled[, intersect(colnames(dat$Y), colnames(pooled))])
    X_pooled <- as.matrix(pooled[, intersect(colnames(dat$X), colnames(pooled))])

    Y_pooled <- 1 * (Y_pooled > 0)

    res$Y_pooled <- Y_pooled
    res$X_pooled <- X_pooled

    if (save_data) {
        saveRDS(Y_pooled, file = file.path(dirs$mod_dat, "Y_pooled.rds"))
        saveRDS(X_pooled, file = file.path(dirs$mod_dat, "X_pooled.rds"))
        saveRDS(res, file = file.path(dirs$mod_dat, "dat.rds"))
    }
    return(res)
            
}
