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
                           totsamp = NA,
                           trans = NA,
                           thn = NA,
                           samps = NA,
                           nchains = NA,
                           nfolds = NA), 
                      class = "varlist")            

    vars$study <- study
    vars$fit <- fit
    vars$totsamp <- sampling$totsamp
    vars$trans <- sampling$trans
    vars$thn <- sampling$thn
    vars$samps <- (sampling$totsamp-sampling$trans)/sampling$thn
    vars$nchains <- sampling$nchains
    vars$nfolds <- sampling$nfolds

    vars$yvars <- c("Clo",
                    "En", 
                    "Be", 
                    "Cap", 
                    "Cau")

    vars$xvars <- switch(fit, 
                         "1" = "(Intercept)", 
                         "2" = "(Intercept)", 
                         "3" = c("Population", 
                                 "Genotype",
                                 "Herbivory",
                                 "Plant.area"),
                         "4" = c("Population", 
                                 "Genotype",
                                 "Herbivory",
                                 "Plant.area"),
                         "5" = c("Population",
                                 "Herbivory",
                                 "Plant.area"),
                         "6" = c("Population",
                                 "Genotype",
                                 "Herbivory",
                                 "Plant.area"))
    vars$boolvars <- switch(fit, 
                             "1" = NA, 
                             "2" = NA, 
                             "3" = "Herbivory",
                             "4" = "Herbivory",
                             "5" = "Herbivory",
                             "6" = "Herbivory")
    vars$dumvars <- switch(fit, 
                           "1" = NA, 
                           "2" = NA, 
                           "3" = c("Population", 
                                   "Genotype"),
                           "4" = c("Population", 
                                   "Genotype"),
                           "5" = "Population",
                           "6" = c("Population", 
                                   "Genotype"))
                                       
    vars$pivars <- switch(fit, 
                          "1" = "Plant", 
                          "2" = "Plant", 
                          "3" = "Plant", 
                          "4" = "Plant", 
                          "5" = "Plant",
                          "6" = "Plant")
                          
    vars$partition <- "Plant"
    
    vars$random <- switch(fit, 
                          "1" = FALSE, 
                          "2" = TRUE, 
                          "3" = FALSE, 
                          "4" = TRUE,
                          "5" = TRUE,
                          "6" = TRUE)

    vars$covDepXvars <- NULL
    if (fit == "6") {    
        vars$covDepXvars <- "Genotype"
        vars$covDepLevel <- 1
    }
    
    return(vars)
}
