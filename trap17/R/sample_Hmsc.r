#' Fit an Hmsc model
#'
#' Fit an Hmsc model with defined variables and save results to defined directories
#' @param param Descr

sample_Hmsc <- function(dat = dat,
                        vars = vars,
                        dirs = dirs,
                        return_ps = TRUE) {

    m1 <- trap17:::prepare_m(dat = dat,
                             vars = vars)

    foldname <- create_name(study = vars$study,
                            totsamp = vars$totsamp,
                            nfolds = vars$nfolds, 
                            type = "fold")
    output_dir <- file.path(dirs$fits, foldname)
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }
    tblname <- paste("model", vars$fit, sep = "_")
    tblname <- paste0(tblname, ".txt")
    vars2 <- vars
    vars2$yvars <- colnames(m1$Y)
    vars2$xvars <- colnames(m1$X)
    vars2$studyDesign <- colnames(m1$studyDesign)
    vars2$random <- m1$ranLevels
    write.table(unlist(vars2), 
                file = file.path(output_dir, tblname), 
                quote = FALSE, 
                sep = ":", 
                col.names = FALSE)

    post_align <- TRUE
    if (!is.null(vars$covDepXvars)) {
        post_align <- FALSE
    }

    print(paste("Model fitted for", vars$study))
    print(paste("Model fit type", vars$fit))
    ps <- NULL
    ps <- Hmsc:::sampleMcmc(m1, 
                            samples = vars$samps,
                            transient = vars$trans,
                            thin = vars$thn,
                            nChains = 2,
                            nParallel = 2,
                            alignPost = post_align)

    filename <- paste("ps", vars$fit, sep = "_")
    filename <- paste0(filename, ".rds")
    saveRDS(ps, file = file.path(output_dir, filename)) 
    if (return_ps) {
        return(ps)
    }
}
