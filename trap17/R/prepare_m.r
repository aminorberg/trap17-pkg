#' Prepare the Hmsc model object
#'
#' Prepare the Hmsc model object for model fitting
#' @param vars Variables, created with set_vars
#' @param dat Data, i.e. trapdata

prepare_m <- function(vars,
                      dat) 
{
    Y <- dat$Y_pooled
    Xorig <- dat$X_pooled
    if (all(vars$xvars == "(Intercept)")) {
        X <- matrix(1, ncol = 1, nrow = nrow(dat$X_pooled))
    } else {
        X <- dat$X_pooled[,vars$xvars]
        Xorig_sel <- X

        if (any(!is.na(vars$dumvars))) {
            dums <- trap17:::dummify(dat = X,
                                     dumvars = vars$dumvars)    
            Xtmp <- data.frame(X[ , !(colnames(X) %in% vars$dumvars)])
            colnames(Xtmp) <- setdiff(colnames(X), vars$dumvars)
            X <- data.frame(Xtmp, dums)
        }
        if (any(!is.na(vars$boolvars))) {
            bools <- trap17:::boolify(dat = Xorig_sel,
                            boolvars = vars$boolvars)    
            Xtmp <- data.frame(X[ , !(colnames(X) %in% vars$boolvars)])
            colnames(Xtmp) <- setdiff(colnames(X), vars$boolvars)
            X <- data.frame(Xtmp, bools)
        }
        if (!is.null(dim(X)) & any(is.na(X))) {
            for (i in 1:ncol(X)) {
                X[which(is.na(X[, i])), i] <- mean(X[, i], na.rm = TRUE)    
            }
        }
        if (is.null(dim(X)) & all(is.na(X))) {
            X <- data.frame(matrix(0, ncol = 1, nrow = nrow(Y)))
            names(X) <- "foo"
        }
    }
    xForm <- switch(vars$fit, 
                    "1" = formula(~1), 
                    "2" = formula(~1),
                    "3" = formula(~.),
                    "4" = formula(~.),
                    "5" = formula(~.),
                    "6" = formula(~.),
                    "7" = formula(~.))

    rLs <- NULL
    if (all(is.na(vars$pivars))) {
        studyDesign <- NULL
    } else {
        studyDesign <- data.frame(Xorig[, vars$pivars])
        for (i in 1:ncol(studyDesign)) {
            studyDesign[, i] <- as.factor(studyDesign[, i])
        }
        names(studyDesign) <- vars$pivars
        if (vars$random) {
            rLs <- list()
            for (i in 1:ncol(studyDesign)) {
                rLs[[i]] <- Hmsc:::HmscRandomLevel(units = studyDesign[,i])
            }
            names(rLs) <- colnames(studyDesign)
            if (!is.na(vars$spat)) {
                spatial_level <- which(colnames(studyDesign) == vars$spat)
                rLs[[spatial_level]] <- Hmsc:::HmscRandomLevel(sData = dat$spat)
            }
        }
        if (!is.null(vars$covDepXvars)) {
            xData <- data.frame(rep(1, nrow(X)), X[, grep("Genotype", colnames(X))])
            colnames(xData) <- c("(Intercept)", colnames(X[, grep("Genotype", colnames(X))]))
            rLs[[vars$covDepLevel]] <- Hmsc:::HmscRandomLevel(xData = xData)        
        }        
    }
    X <- as.data.frame(X)

    m1 <- Hmsc:::Hmsc(Y = as.matrix(Y), 
                      XData = X,
                      XFormula = xForm,
                      studyDesign = studyDesign,
                      ranLevels = rLs,
                      distr = "probit",
                      XScale = TRUE,
                      XSelect = NULL,
                      phyloTree = NULL,
                      C = NULL)
    return(m1)
}
