#' Compute associations (modified)
#'
#' Modified function for computing association matrices when using covariate dependent latent variables.
#' @param object A fitted Hmsc model object.
#' @param ... See ?Hmsc::Compute associations for the rest
#' @export

computeAssociations_modified <- function (hM, start = 1, thin = 1) 

{
    OmegaCor = vector("list", hM$nr)
    postList = poolMcmcChains(hM$postList, start = start, thin = thin)
    getOmegaCor <- function(a, r = r) {
        tmp <- list()
        for (i in 1:dim(a$Lambda[[r]])[3]) {
            tmp[[i]] <- cov2cor(crossprod(a$Lambda[[r]][,,i]))
        }
        res <- simplify2array(tmp)
        return(res)
    } 

    for (r in seq_len(hM$nr)) {
        OmegaCor1 = lapply(postList, getOmegaCor, r = r)
        mOmegaCor1 = apply(abind::abind(simplify2array(OmegaCor1)), 
                           c(1, 2, 3), 
                           mean)
        OmegaCor2 = lapply(OmegaCor1, function(a) return(a > 0))
        support1 = apply(abind::abind(simplify2array(OmegaCor2)), 
                         c(1, 2, 3), 
                         mean)
        dimnames(mOmegaCor1) <- list(hM$spNames, hM$spNames, colnames(hM$rL[[r]]$x))
        dimnames(support1) <- list(hM$spNames, hM$spNames, colnames(hM$rL[[r]]$x))
        tmp = list()
        tmp$mean = mOmegaCor1
        tmp$support = support1
        OmegaCor[[r]] = tmp
    }
    return(OmegaCor)
}

