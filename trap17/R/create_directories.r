#' Create directories
#'
#' Create the folders/directories in which model fits and results will be saved.
#' @param foldname A 'trapfilename' object created with create_name, specifying the name for the study folder.
#' @param dirs A 'dirlist' object created with set_dirs.
#' @export

create_directories <- function(dirs,
                               foldname)
{
    dir.create(file.path(dirs$fits, foldname))
    dir.create(file.path(dirs$fits, foldname, "figs"))

}
