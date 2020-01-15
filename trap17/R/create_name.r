#' Create name
#'
#' Create a name for a folder or a cross-validation results file 
#' @param study Name of the study
#' @param totsamp Total number of MCMC samples
#' @param nfolds Number of folds for cross-validation
#' @param obj_type Object type: either "fold", "cv" or "eval_cv"
#' @export

create_name <- function(study = "trap17", 
                        totsamp = NULL,
                        nfolds = NULL, 
                        type = c("fold", "cv", "eval_cv"))
{
    filebody <- paste0(study, 
                       "_totsamp",
                       totsamp)
    res <- switch(type,
                  fold = paste0(study, "_totsamp", totsamp),
                  cv = paste0("cv_preds_", "nfolds", nfolds),
                  eval_cv = paste0("eval_cv_", "nfolds", nfolds))
    class(res) <- "trapfilename"
    return(res)
}
