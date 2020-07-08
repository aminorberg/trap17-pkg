#' Perform cross-validation
#'
#' Perform cross-validation
#' @param param Descr

do_cv <- function(ps,
                  dat,
                  dirs,
                  vars,
                  higher_eval_levels = TRUE, 
                  save_cv = TRUE,
                  start_iter,
                  expect = "both")
{

    eval_cv <- NA
    pred_corrs <- NA
    cv_preds <- NA
    cv_preds_realz <- NA

    if (!is.character(expect)) {
        stop("The parameter 'expected' has to be a character, 'TRUE'/'FALSE'/'both'")
    }
    expectations <- switch(expect,
                           "TRUE" = c(TRUE, FALSE),
                           "FALSE" = c(FALSE, TRUE),
                           "both" = c(TRUE, TRUE))

    res <- structure(list(evaluation_cv = NULL, 
                          higher_evaluation_cv = NULL, 
                          cv_predictions = NULL, 
                          cv_predictions_realisations = NULL),
                     class = "cvresults")
    
    set.seed(7)
    cv_partition <- Hmsc:::createPartition(hM = ps, 
                                           nfolds = vars$sampling$nfolds, 
                                           column = vars$partition)

    if (!is.null(vars$covDepXvars)) {    
        if (expectations[1]) {
            cv_preds <- trap17:::computePredictedValues_modified(hM = ps, 
                                                                 partition = cv_partition, 
                                                                 expected = TRUE,
                                                                 start = start_iter,
                                                                 alignPost = FALSE)
        }
        if (expectations[2]) {
            cv_preds_realz <- trap17:::computePredictedValues_modified(hM = ps, 
                                                                       partition = cv_partition, 
                                                                       expected = FALSE,
                                                                       start = start_iter,
                                                                       alignPost = FALSE)
        }
    } else {
        if (expectations[1]) {
            cv_preds <- Hmsc:::computePredictedValues(hM = ps, 
                                                      partition = cv_partition, 
                                                      expected = TRUE,
                                                      start = start_iter,
                                                      alignPost = TRUE)
        }
        if (expectations[2]) {
            cv_preds_realz <- Hmsc:::computePredictedValues(hM = ps, 
                                                            partition = cv_partition, 
                                                            expected = FALSE,
                                                            start = start_iter,
                                                            alignPost = TRUE)
        }
    }

    if (save_cv) {
        foldname <- create_name(study = vars$study,
                            totsamp = vars$sampling$totsamp,
                            nfolds = vars$sampling$nfolds, 
                            type = "fold")

        output_dir <- file.path(dirs$fits, foldname)
        if (!dir.exists(output_dir)) {
            dir.create(output_dir)
        }

        cv_filename <- create_name(totsamp = vars$sampling$totsamp,
                                   nfolds = vars$sampling$nfolds, 
                                   type = "cv")
        cv_filename <- paste(cv_filename, "ps", vars$fit, sep = "_")
        if (expectations[1]) {
            saveRDS(cv_preds, 
                    file = file.path(output_dir, paste0(cv_filename, ".rds")))
        }
        if (expectations[2]) {
            saveRDS(cv_preds_realz, 
                    file = file.path(output_dir, paste0(cv_filename, "_realz.rds")))
        }
    }
    if (expectations[1]) {
        eval_cv <- Hmsc:::evaluateModelFit(hM = ps, predY = cv_preds)
    }

    if (higher_eval_levels) {
        if (expectations[2]) {
            cv_preds_for_high <- cv_preds_realz
            print("Higher level performance is calculated for realisations")
        } else {
            cv_preds_for_high <- cv_preds    
            print("Higher level performance is calculated for expected values")
        }    
        pred_corrs <- trap17:::higher_cors(dat = dat,
                                           preds = cv_preds_for_high)
    }

    if (save_cv) {
        cv_eval_filename <- create_name(nfolds = vars$sampling$nfolds, 
                                        type = "eval_cv")
        cv_eval_filename <- paste(cv_eval_filename, "ps", vars$fit, sep = "_")
        saveRDS(eval_cv, 
                file = file.path(output_dir, paste0(cv_eval_filename, ".rds")))
        if (!is.null(higher_eval_levels)) {
            saveRDS(pred_corrs, 
                    file = file.path(output_dir, paste0("higher_",
                                                               cv_eval_filename,
                                                               ".rds")))
        }
    }
    res$cv_predictions <- cv_preds
    res$cv_predictions_realisations <- cv_preds_realz
    res$higher_evaluation_cv <- pred_corrs
    res$evaluation_cv <- eval_cv
    
    return(res)

}
