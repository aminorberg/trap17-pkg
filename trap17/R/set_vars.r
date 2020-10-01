#' Set variables for the model
#'
#' Load objects from a directory: model fit, predictions, cross-validation results, ...
#' @param study Name of study
#' @param fit Model fit
#' @param sampling a "mcmcsettings" object defining the sampling settings
#' @export

set_vars <- function(study = "trap17",
                     fit,
                     sampling) {

    if (class(sampling) != "mcmcsettings") {
        stop((print("parameter 'sampling' is not a  'mcmcsettings' object")))
    }

    vars <- structure(list(study = NA,
                           fit = NA,
                           xvars = NA,
                           pivars = NA,
                           studyDesign = NA,
                           partition = NA,
                           random = FALSE,
                           spat = NA,
                           boolvars = NA,
                           dumvars = NA,
                           covDepXvars = NA,
                           covDepLevel = NA,
                           sampling = NA), 
                      class = "varlist")            

    vars$study <- study
    vars$fit <- fit
    vars$sampling <- sampling
    vars$sampling$samps <- (sampling$totsamp-sampling$trans)/sampling$thn
    vars$sampling$mod_rl_priors <- sampling$mod_rl_priors

    vars$yvars <- c("Clo",
                    "Be", 
                    "Cap", 
                    "Cau",
                    "En")

    vars$xvars <- switch(fit, 
                         "1" = c("Population",
                                 "Herbivory",
                                 "Plant.area"),
                         "2" = c("Population", 
                                 "Genotype",
                                 "Herbivory",
                                 "Plant.area"),
                         "3" = c("Population",
                                 "Genotype",
                                 "Herbivory",
                                 "Plant.area"))
    vars$boolvars <- "Herbivory"
    vars$dumvars <- switch(fit, 
                           "1" = "Population",
                           "2" = c("Population", 
                                   "Genotype"),
                           "3" = c("Population", 
                                   "Genotype"))
                                       
    vars$pivars <- "Plant"
    vars$partition <- "Plant"
    
    vars$random <- TRUE

    vars$covDepXvars <- NULL
    if (fit == "3") {    
        vars$covDepXvars <- "Genotype"
        vars$covDepLevel <- 1
    }
    
    return(vars)
}
