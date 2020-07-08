#' Do cross-validation
#'
#' Loads existing model objects and does cross-validation. If do_cv = TRUE, returns the performance measures.
#' @param dat Data
#' @param dirs Directories
#' @param variants Selection of model variants
#' @param sampling Sampling settings
#' @param saveCVs default TRUE
#' @export

just_cv <- function(dat, 
                    dirs,
                    variants = "ALL",
                    sampling,
                    start_iter,
                    expect = "both",
                    saveCVs = TRUE) 
{

    if (any(variants == "ALL")) {
        fitss <- as.character(1:3)
    } else {
        fitss <- as.character(variants)
    }
    res <- list()

    for (f in 1:length(fitss)) {

        vars <- trap17:::set_vars(study = "trap17",
                                  fit = fitss[f],
                                  sampling = sampling)
        foldname <- create_name(study = vars$study,
                                totsamp = vars$sampling$totsamp,
                                nfolds = vars$sampling$nfolds, 
                                type = "fold")
        output_dir <- file.path(dirs$fits, foldname)
        filename <- paste("ps", vars$fit, sep = "_")
        filename <- paste0(filename, ".rds")
        ps <- NULL
        ps <- readRDS(file = file.path(output_dir, filename)) 
        res[[f]] <- trap17:::do_cv(ps = ps,
                                   dat = dat,
                                   dirs = dirs,
                                   vars = vars,
                                   start_iter = start_iter,
                                   expect = expect,
                                   save_cv = TRUE,
                                   higher_eval_levels = TRUE)               
    }                     
    names(res) <- paste0("ps", fitss)
    return(res)
}
