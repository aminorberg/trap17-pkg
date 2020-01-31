#' Fit models and do cross-validation
#'
#' Fit models and optionally does cross-validation. If do_cv = TRUE, returns the performance measures.
#' @param dat Data
#' @param dirs Directories
#' @param variants Selection of model variants
#' @param sampling Sampling settings
#' @param saveCVs default TRUE
#' @export

model_and_cv <- function(dat, 
                         dirs,
                         variants = "ALL",
                         sampling,
                         do_cv = TRUE,
                         saveCVs = TRUE) 
{

    if (any(variants == "ALL")) {
        fitss <- as.character(1:7)
    } else {
        fitss <- as.character(variants)
    }
    res <- list()

    for (f in 1:length(fitss)) {

        vars <- trap17:::set_vars(study = "trap17",
                                  fit = fitss[f],
                                  sampling = sampling)

        ps <- trap17:::sample_Hmsc(dat = dat,
                                   vars = vars,
                                   dirs = dirs,
                                   return_ps = TRUE)
        if (do_cv) {
            res[[f]] <- trap17:::do_cv(ps = ps,
                                       dat = dat,
                                       dirs = dirs,
                                       vars = vars,
                                       save_cv = TRUE,
                                       higher_eval_levels = TRUE)               
        } 
    }                     
    if (do_cv) {
        names(res) <- paste0("ps", fitss)
        if (saveCVs) {
            foldname <- create_name(study = vars$study,
                                    totsamp = vars$totsamp,
                                    nfolds = vars$nfolds, 
                                    type = "fold")
            output_dir <- file.path(dirs$fits, foldname)
            filename <- "cv_res_all.rds"
            saveRDS(res, file = file.path(output_dir, filename))
        }
        return(res)
    }
}
