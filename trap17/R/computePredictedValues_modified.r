#' Compute predicted values (modified)
#'
#' Modified function for computing predicted values, using the predictHmsc_modified function for predictions
#' @param hM Hmsc object
#' @param partition Partition for the cross validation
#' @param alignPost boolean

#partition.sp = NULL
#start = 1
#thin = 1
#Yc = NULL
#mcmcStep = 1
#initPar = NULL 
#nParallel = 1
#nChains = length(hM$postList)
#updater = list()
#verbose = hM$verbose

computePredictedValues_modified <- function (hM, partition = NULL, partition.sp = NULL, start = 1, 
                                             thin = 1, Yc = NULL, mcmcStep = 1, expected = TRUE, initPar = NULL, 
                                             nParallel = 1, nChains = length(hM$postList), updater = list(), 
                                             verbose = hM$verbose, alignPost = TRUE) 
{


    if (is.null(partition)) {
        postList = poolMcmcChains(hM$postList, start = start, 
            thin = thin)
        pred = predict(hM, post = postList, Yc = Yc, mcmcStep = 1, 
            expected = expected)
        predArray = abind(pred, along = 3)
    } else {
        if (length(partition) != hM$ny) {
            stop("HMSC.computePredictedValues: partition parameter must be a vector of length ny")
        }
        nfolds = length(unique(partition))
        postN = Reduce(sum, lapply(hM$postList, length))
        predArray = array(NA, c(hM$ny, hM$ns, postN))

        for (k in 1:nfolds) {
            print(sprintf("Cross-validation, fold %d out of %d", k, nfolds))
            train = (partition != k)
            val = (partition == k)
            dfPi = as.data.frame(matrix(NA, sum(train), hM$nr))
            colnames(dfPi) = hM$rLNames
            for (r in seq_len(hM$nr)) {
                dfPi[, r] = factor(hM$dfPi[train, r])
            }
            switch(class(hM$X)[1L], matrix = {
                XTrain = hM$X[train, , drop = FALSE]
                XVal = hM$X[val, , drop = FALSE]
            }, list = {
                XTrain = lapply(hM$X, function(a) a[train, , 
                  drop = FALSE])
                XVal = lapply(hM$X, function(a) a[val, , drop = FALSE])
            })
            if (hM$ncRRR > 0) {
                XRRRTrain = hM$XRRR[train, , drop = FALSE]
                XRRRVal = hM$XRRR[val, , drop = FALSE]
            } else {
                XRRRTrain = NULL
                XRRRVal = NULL
            }
            hM1 = Hmsc(Y = hM$Y[train, , drop = FALSE], X = XTrain, 
                XRRR = XRRRTrain, ncRRR = hM$ncRRR, XSelect = hM$XSelect, 
                distr = hM$distr, studyDesign = dfPi, Tr = hM$Tr, 
                C = hM$C, ranLevels = hM$rL)
            setPriors(hM1, V0 = hM$V0, f0 = hM$f0, mGamma = hM$mGamma, 
                UGamma = hM$UGamma, aSigma = hM$aSigma, bSigma = hM$bSigma, 
                nu = hM$nu, a1 = hM$a1, b1 = hM$b1, a2 = hM$a2, 
                b2 = hM$b2, rhopw = hM$rhowp)
            hM1$YScalePar = hM$YScalePar
            hM1$YScaled = (hM1$Y - matrix(hM1$YScalePar[1, ], 
                hM1$ny, hM1$ns, byrow = TRUE))/matrix(hM1$YScalePar[2, 
                ], hM1$ny, hM1$ns, byrow = TRUE)
            hM1$XInterceptInd = hM$XInterceptInd
            hM1$XScalePar = hM$XScalePar
            switch(class(hM$X)[1L], matrix = {
                hM1$XScaled = (hM1$X - matrix(hM1$XScalePar[1, 
                  ], hM1$ny, hM1$ncNRRR, byrow = TRUE))/matrix(hM1$XScalePar[2, 
                  ], hM1$ny, hM1$ncNRRR, byrow = TRUE)
            }, list = {
                hM1$XScaled = list()
                for (zz in seq_len(length(hM1$X))) {
                  hM1$XScaled[[zz]] = (hM1$X[[zz]] - matrix(hM1$XScalePar[1, 
                    ], hM1$ny, hM1$ncNRRR, byrow = TRUE))/matrix(hM1$XScalePar[2, 
                    ], hM1$ny, hM1$ncNRRR, byrow = TRUE)
                }
            })
            if (hM1$ncRRR > 0) {
                hM1$XRRRScalePar = hM$XRRRScalePar
                hM1$XRRRScaled = (hM1$XRRR - matrix(hM1$XRRRScalePar[1, 
                  ], hM1$ny, hM1$ncORRR, byrow = TRUE))/matrix(hM1$XRRRScalePar[2, 
                  ], hM1$ny, hM1$ncORRR, byrow = TRUE)
            }
            hM1$TrInterceptInd = hM$TrInterceptInd
            hM1$TrScalePar = hM$TrScalePar
            hM1$TrScaled = (hM1$Tr - matrix(hM1$TrScalePar[1, 
                ], hM1$ns, hM1$nt, byrow = TRUE))/matrix(hM1$TrScalePar[2, 
                ], hM1$ns, hM1$nt, byrow = TRUE)
            hM1 = sampleMcmc(hM1, samples = hM$samples, thin = hM$thin, 
                transient = hM$transient, adaptNf = hM$adaptNf, 
                initPar = initPar, nChains = nChains, nParallel = nParallel, 
                updater = updater, verbose = verbose, alignPost = alignPost)
            postList = poolMcmcChains(hM1$postList, start = start)
            dfPi = as.data.frame(matrix(NA, sum(val), hM$nr))
            colnames(dfPi) = hM$rLNames
            for (r in seq_len(hM$nr)) {
                dfPi[, r] = factor(hM$dfPi[val, r])
            }

            if (is.null(partition.sp)) {
                pred1 = trap17:::predictHmsc_modified(hM1, post = postList, X = XVal, 
                                                      XRRR = XRRRVal, studyDesign = dfPi, 
                                                      Yc = Yc[val, , drop = FALSE], 
                                                      mcmcStep = mcmcStep, expected = expected)
#               pred1 = predict(hM1, post = postList, X = XVal, 
#                               XRRR = XRRRVal, studyDesign = dfPi, 
#                               Yc = Yc[val, , drop = FALSE], 
#                               mcmcStep = mcmcStep, expected = expected)
                pred1Array = abind::abind(pred1, along = 3)
            } else {
                pred1Array = array(dim = c(sum(val), hM$ns, postN))
                nfolds.sp = length(unique(partition.sp))
                for (i in 1:nfolds.sp) {
                  train.sp = (partition.sp != i)
                  val.sp = (partition.sp == i)
                  Yc = matrix(NA, nrow = sum(val), ncol = hM$ns)
                  Yc[, train.sp] = hM$Y[val, train.sp, drop = FALSE]
                  pred2 = predict(hM1, post = postList, X = XVal, 
                    XRRR = XRRRVal, studyDesign = dfPi, Yc = Yc, 
                    mcmcStep = mcmcStep, expected = expected)
                  pred2Array = abind(pred2, along = 3)
                  pred1Array[, val.sp, ] = pred2Array[, val.sp, 
                    ]
                }
            }
            predArray[val, , ] = pred1Array
        }
    }
    return(predArray)
}
