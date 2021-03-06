#' Sampling settings
#'
#' Set settings for MCMC sampling with Hmsc
#' @param totsamp Number of MCMC samples
#' @param trans Number of transition MCMC samples
#' @param thn Thinning of MCMC samples
#' @param nfolds Number of folds for cross-validation
#' @return List of settings to be passed on for MCMC sampling
#' @export

sampling_settings <- function (totsamp = 100,
                               trans = 50,
                               thn = 1,
                               nfolds = 2,
                               nchains = 2,
                               mod_rl_priors = NA)
{
    res <- structure(list(totsamp = NULL, 
                          trans = NULL, 
                          thn = NULL, 
                          nchains = NULL, 
                          nfolds = NULL), 
                     class = "mcmcsettings")
    res$totsamp <- totsamp
    res$trans <- trans
    res$thn <- thn
    res$nchains <- nchains
    res$nfolds <- nfolds
    res$mod_rl_priors <- mod_rl_priors
    
    return(res)
}

