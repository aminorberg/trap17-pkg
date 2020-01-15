#' Select data
#'
#' Select data
#' @param vars Descr
#' @param dirs Descr
#' @export

sel_data <- function(vars,
                     dirs) {

    res <- structure(list(Y = NULL, X = NULL, PI = NULL, studyDesign = NULL), 
                     class = "trapdata")

    if (vars$study == "trap") {
        dat <- readRDS(file = file.path(dirs$mod_dat, "dat.rds"))
    }
    
    spat <- NA
    if (!is.na(vars$spat)) {
        spat <- readRDS(file = file.path(dirs$mod_dat, "spat.rds"))
    }
    
    dfr <- data.frame(dat$Y_pooled, dat$X_pooled, dat$PI_pooled)
    Y <- data.frame(dfr[,vars$yvars])
    
    X <- NA
    if (all(!is.na(vars$xvars))) {
        X <- data.frame(dfr[,vars$xvars])
        names(X) <- vars$xvars
    }
    PI <- NA
    if (all(!is.na(vars$pivars))) {
        PI <- data.frame(dfr[,vars$pivars])
        names(PI) <- vars$pivars
    }
    studyDesign <- data.frame(dfr[,vars$studyDesign])
    
    res$Y <- Y
    res$X <- X
    res$PI <- PI
    res$spat <- spat
    res$studyDesign <- studyDesign

    return(res)
}
