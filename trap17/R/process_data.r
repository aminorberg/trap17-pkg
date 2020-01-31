#' Process trap17 data
#'
#' Modify variable names and pool the data
#' @param filename Original data file
#' @param dirs A trap17 class "dirlist" object
#' @param save_data Save the modified data
#' @param return_data Return the pooled data
#' @export


process_data <- function(filename,
                         dirs,
                         save_data = TRUE, 
                         return_data = TRUE)
{

    full_data <- data_prep(filename = filename,
                           save_data = save_data,
                           dirs = dirs)

    dat_pooled <- pool(dat = full_data,
                       save_data = save_data,
                       dirs = dirs)

    if (return_data) {
        full_data$Y_pooled <- dat_pooled$Y_pooled
        full_data$X_pooled <- dat_pooled$X_pooled

        return(full_data)
    }
    
}
