#' Data preparation
#'
#' Prepare the original data: unify names and define aggregation formulas
#' @param filename Original data file
#' @param dirs A trap17 class "dirlist" object

data_prep <- function(filename = "TRAP17.csv",
                      dirs = dirs) {

    dat <- structure(list(Y = NULL,
                          Y_pooled = NULL,
                          X = NULL,
                          X_pooled = NULL,
                          PI = NULL, 
                          PI_pooled = NULL, 
                          yvars = NULL, 
                          xvars = NULL,
                          pivars = NULL,
                          pool_sum_formula = NULL,
                          pool_min_formula = NULL,
                          pool_max_formula = NULL,
                          grouping_var = NULL), 
                     class = "trapdata")

    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    dat_tmp <- read.csv2(file.path(datapath)) #check.names = FALSE
    dat_tmp <- data.frame(dat_tmp)
    
    names(dat_tmp)[which(names(dat_tmp) == "Clostero")] <- "Clo"
    names(dat_tmp)[which(names(dat_tmp) == "Polero")] <- "En"
    names(dat_tmp)[which(names(dat_tmp) == "Partiti")] <- "Be"
    names(dat_tmp)[which(names(dat_tmp) == "PILV")] <- "Pl"
    names(dat_tmp)[which(names(dat_tmp) == "CAMV")] <- "Ca"

    
    dat$yvars <- c("Clo","En", "Be", "Pl", "Ca")
    dat$xvars <- c("Genotype", "Population", "Plant.area", "Herbivory")
    dat$pivars <- c("sampleID", "Timepoint", "Plant")
    dat$pivars_pooled <- "Plant"
    dat$pool_sum_formula <- formula(cbind(Clo, 
                                          En, 
                                          Be, 
                                          Pl, 
                                          Ca, 
                                          Herbivory, 
                                          Timepoint) ~ Plant)
    dat$pool_min_formula <- formula(cbind(Genotype, Population) ~ Plant)
    dat$pool_max_formula <- formula(Plant.area ~ Plant)
    dat$grouping_var <- "Plant"

    Ydat <- dat_tmp[, dat$yvars]
    Ydat <- apply(Ydat, 2, as.numeric)
    if (!all(apply(Ydat, 2, is.numeric))) {
        stop("Y is not numeric")
    }

    Xdat <- dat_tmp[, dat$xvars]
    Xdat <- apply(Xdat, 2, as.numeric)
    Xdat <- data.frame(Xdat)
    if (!all(apply(Xdat, 2, is.numeric))) {
        stop("X is not numeric")
    }
    PIdat <- dat_tmp[, dat$pivars]
    if  ( any(unlist(strsplit(as.character(PIdat[1, "Plant"]), split = "")) == "S") ) {
        PIdat[, "Plant"] <- sub(pattern = "S", 
                                replacement = "", 
                                x = PIdat[, "Plant"])
    }
    PIdat <- apply(PIdat, 2, as.numeric)
    if (!all(apply(PIdat, 2, is.numeric))) {
        stop("PI is not numeric")
    }

    saveRDS(Ydat, file = file.path(dirs$mod_dat, "Y.rds"))
    saveRDS(Xdat, file = file.path(dirs$mod_dat, "X.rds"))
    saveRDS(PIdat, file = file.path(dirs$mod_dat, "PI.rds"))

    dat$Y <- Ydat
    dat$X <- Xdat
    dat$PI <- PIdat
    
    return(dat)

}
