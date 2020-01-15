#' Load objects
#'
#' Load objects from a directory: model fit, predictions, cross-validation results, ...
#' @param path Path to object
#' @param study Name of study
#' @param obj_type Object type: either "fold", "cv" or "eval_cv"
#' @export

load_objects_from_dir <- function(path, 
                                  study,
                                  obj_type) {

    res <- structure(list(), class = "objlist")            
    if (obj_type == "ps") {
            obj_type <- paste0("^", obj_type)
    } 
    fold <- list.files(path = path, 
                       pattern = study, 
                       full.names = FALSE, 
                       ignore.case = FALSE)
            
    if (length(fold) > 1) {
        stop("Too many folders with tho study name behind the specified path")
    }
    obj <- NULL
    obj <- list.files(path = file.path(path, fold), 
                      pattern = obj_type, 
                      full.names = FALSE, 
                      ignore.case = FALSE)            

    for (k in 1:length(obj)) {
        foc_file <- file.path(path, fold, obj[k])
        
        if (file.exists(foc_file)) {
            loaded <- readRDS(file = foc_file)    
            res <- append(res, list(loaded))
            names(res)[length(res)] <- paste(fold, obj[k], sep = "_")            
        }
    }
    return(res)
}
