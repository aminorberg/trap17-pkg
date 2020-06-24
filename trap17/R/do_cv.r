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
                  expected = TRUE)
{

    res <- structure(list(eval_cv = NULL, higher_eval_cv = NULL),
                     class = "cvresults")
  
    cv_partition <- Hmsc:::createPartition(hM = ps, 
                                           nfolds = vars$sampling$nfolds, 
                                           column = vars$partition)

    if (!is.null(vars$covDepXvars)) {
        cv_preds <- trap17:::computePredictedValues_modified(hM = ps, 
                                                             partition = cv_partition, 
                                                             expected = expected,
                                                             start = start_iter,
                                                             alignPost = FALSE)
    } else {
        cv_preds <- Hmsc:::computePredictedValues(hM = ps, 
                                                  partition = cv_partition, 
                                                  expected = expected,
                                                  start = start_iter,
                                                  alignPost = TRUE)
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
        saveRDS(cv_preds, file = file.path(output_dir, paste0(cv_filename, ".rds")))
    }
    eval_cv <- Hmsc:::evaluateModelFit(hM = ps, predY = cv_preds)

    if (higher_eval_levels) {
        form <- formula(cbind(Clo, 
                              Be, 
                              Cap, 
                              Cau, 
                              En) ~ Population + Genotype)

        pooled_occ_true <- as.matrix(aggregate(formula = form, 
                                               data = cbind(dat$X_pooled, dat$Y_pooled), 
                                               FUN = sum))[, -c(1:2)]
        corrs <- matrix(NA, 
                        nrow = dim(cv_preds)[3],
                        ncol = ncol(ps$Y))
        for (j in 1:dim(cv_preds)[3]) {
            tmp <- cv_preds[, , j]
            colnames(tmp) <- colnames(ps$Y)
            tmp <- cbind(tmp, dat$X_pooled)
            preds <- as.matrix(aggregate(formula = form, 
                                         data = tmp, 
                                         FUN = sum))[, -c(1:2)]
            corrs[j,] <- diag(apply(preds, 2, cor, pooled_occ_true, method = "spearman"))
        }
    }

    if (save_cv) {
        cv_eval_filename <- create_name(nfolds = vars$nfolds, 
                                        type = "eval_cv")
        cv_eval_filename <- paste(cv_eval_filename, "ps", vars$fit, sep = "_")
        saveRDS(eval_cv, 
                file = file.path(output_dir, paste0(cv_eval_filename, ".rds")))
        if (!is.null(higher_eval_levels)) {
            saveRDS(corrs, 
                    file = file.path(output_dir, paste0("higher_",
                                                               cv_eval_filename,
                                                               ".rds")))
        }
    }
    res$cv_predictions <- cv_preds
    res$higher_eval_cv <- corrs
    res$eval_cv <- eval_cv
    
    return(res)

}
