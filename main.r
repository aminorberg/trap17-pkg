# TRAP17 ANALYSIS 2020

# 1 DATA: importation and raw data figures
rm(list = ls(all = TRUE)) ; gc()
working_dir <- "/Users/anorberg/Documents/Zurich/UZH/TRAP/pkg/trap17-pkg"
setwd(working_dir)
library("trap17")
dirs <- set_dirs(working_dir = working_dir)
#saveRDS(dirs, file = file.path(working_dir, "dirs.rds"))

# 1.1a process...
#dat <- process_data(filename = "TRAP17.csv",
#                  dirs = dirs,
#                  save_data = TRUE,
#                  return_data = TRUE)

# 1.1b ...or load the processed data
dat <- trapdata
str(dat)

# 2 MODEL FITTING

# 2.1 sampling settings
sampling <- sampling_settings(totsamp = 1000,
                              trans = 500,
                              thn = 1,
                              nchains = 1,
                              nfolds = 3)
foldname <- create_name(study = "trap17",
                        totsamp = sampling$totsamp,
                        nfolds = sampling$nfolds, 
                        type = "fold")
create_directories(foldname = foldname, dirs = dirs) # muokkaa yleisemmäksi
saveRDS(sampling, 
        file = file.path(dirs$fits, foldname, "sampling.rds"))

# 2.2 model fitting and cross-validation
evals <- model_and_cv(dat = dat,
                      dirs = dirs,
                      variants = 1:7,
                      sampling = sampling,
                      saveCVs = TRUE)

str(evals)
saveRDS(evals, 
        file = file.path(dirs$fits, foldname, "evals.rds"))

# 3 RESULTS
#rm(list = setdiff(ls(), c("working_dir", "dirs", "foldname", "sampling", "dat"))); gc()
#setwd(working_dir)
#library("trap17")
foc_study <- "trap17_totsamp1e"
foldname <- "trap17_totsamp1e+05"

evals <- readRDS(file = file.path(file.path(dirs$fits, foldname), "evals.rds"))

# 3.1 cv-based R2s
tjurs <- lapply(lapply(evals, '[[', 1), '[[', 3)
cors <- lapply(lapply(evals, '[[', 2), colMeans, na.rm = TRUE)

lapply(tjurs, mean)
lapply(cors, mean)
ps_sel <- c(3:4, 7)
lapply(tjurs, mean)[ps_sel]
lapply(cors, mean)[ps_sel]
lapply(tjurs, sd)[ps_sel]
lapply(cors, sd)[ps_sel]

# 3.2 Co-occurrence combinations
orig_combs <- co_occ_combs(partition = dat$X_pooled[,c("Genotype","Population")],
                           Y_arr = array(dat$Y_pooled, 
                                         dim = c(nrow(dat$Y_pooled), 
                                                 ncol(dat$Y_pooled), 
                                                 1),
                                         dimnames = list(1:nrow(dat$Y_pooled),
                                                         colnames(dat$Y_pooled),
                                                         1))) 
saveRDS(orig_combs, 
        file = file.path(dirs$fits, foldname, "orig_combs.rds"))

modelled_combs <- vector(mode = "list", length = length(evals))
for (i in 1:length(evals)) {
    yarr <- evals[[i]]$cv_predictions
    dimnames(yarr) <- list(1:dim(yarr)[1],
                           colnames(dat$Y_pooled),
                           1:dim(yarr)[3])
    modelled_combs[[i]] <- co_occ_combs(partition = dat$X_pooled[,c("Genotype",
                                                                     "Population")],
                                         Y_arr = yarr) 
}
saveRDS(modelled_combs, 
        file = file.path(dirs$fits, foldname, "modelled_combs.rds"))

str(orig_combs)
str(modelled_combs)

all_virus_combs <- as.character(unique(unlist(lapply(lapply(modelled_combs, 
                                                            dimnames), 
                                                            '[', 
                                                            1))))

all_virus_combs <- all_virus_combs[order(sapply(all_virus_combs, nchar), decreasing = TRUE)]
all_virus_combs <- c("Empty", all_virus_combs[-which(all_virus_combs == "Empty")])


prev_ord_pop <- c(3, 2, 4, 1)
prev_ord_genot <- c(2, 4, 1, 3)


colrs <- load_colour_palette()
colmat <- cbind(all_virus_combs, colrs[[2]])
rownames(colmat) <- colmat[,1]


# 3.2.1 plot original co-occurrence combinations
pdf(file = file.path(dirs$fits, 
                     foldname, 
                     "figs",
                     "orig_cooccs.pdf"),
    bg = "transparent", 
    width = 15, 
    height = 5)
    par(family = "serif", mfrow = c(1,4))
    for (g in prev_ord_genot) {
        tmp11 <- c()
        for (p in prev_ord_pop) {
            tmp12 <- orig_combs[,g,p][rownames(colmat)]
            tmp11 <- cbind(tmp11, tmp12)
        }
    if (any(is.na(tmp11))) { 
        tmp11 <- tmp11[-which(is.na(tmp11), arr.ind = TRUE)[,1],]
    }
    colrs1 <- colmat[rownames(tmp11), 2]
    colrs1 <- colrs1[nrow(tmp11):1]
    barplot(tmp11[nrow(tmp11):1,],
            col = colrs1,
            xaxt = "n")
    }
dev.off()

# 3.2.2 predicted co-occurrence combinations
for (i in c(1:length(modelled_combs))) {
    pdf(file = file.path(dirs$fits, 
                         foldname, 
                         "figs",
                         paste0("co_occs_ps", 
                                i,
                                ".pdf")
                        ),
        bg = "transparent", 
        width = 15, 
        height = 5)
    par(family = "serif", mfrow = c(1,4))
    for (g in prev_ord_genot) {
        tmp12 <- NULL
        tmp11 <- c()
        for (p in prev_ord_pop) {
                tmp12 <- modelled_combs[[i]][,g,p][rownames(colmat)]
                tmp11 <- cbind(tmp11, tmp12)
        }
        if (any(is.na(tmp11))) { 
            tmp11 <- tmp11[-which(is.na(tmp11), arr.ind = TRUE)[,1],]
        }
        colrs1 <- colmat[rownames(tmp11), 2]
        colrs1 <- colrs1[nrow(tmp11):1]
        barplot(tmp11[nrow(tmp11):1,],
                col = colrs1,
                xaxt = "n")
                #legend.text = rownames(tmp11)[nrow(tmp11):1])
    }
    dev.off()
}

# 3.2.3 legend
pdf(file = file.path(dirs$fits, 
                     foldname, 
                     "figs", 
                     "cooccurrencebars_legend.pdf"),
    bg = "transparent", 
    width = 5, 
    height = 10)
    par(family = "serif")
    plot(x = 1:10, 1:10, type = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("topleft", legend = colmat[,1], fill = colmat[,2], bty = 'n')
dev.off()

# 3.3 results for the best model variant
pss <- load_objects_from_dir(path = dirs$fits, 
                             study = foc_study,
                             obj_type = "ps")
names(pss)
whichPs <- 3
ps <- pss[[whichPs]]

library(wesanderson)
wesandcols <- wes_palette("Cavalcanti1")[5:1]

# 3.3.1 variance partitioning
if (whichPs > 2 & whichPs != 7) {
    sel <-  c(1, 2, rep(3, 3), rep(4, 3), 5) 
    selnames <- c("Intercept", 
                  "Plant.area", 
                  "Population", 
                  "Genotype", 
                  "Herbivory")

    if (whichPs == 5 | whichPs == 6) {
        selnames <- selnames[-which(selnames == "Genotype")]
        sel <- c(1, 2, rep(3, 3), 4)
    }
    
    vp <- Hmsc::computeVariancePartitioning(ps, 
                                      group = sel, 
                                      groupnames = selnames)
    toPlot <- vp$vals[-1,]
    toPlot <- toPlot[nrow(toPlot):1,]
    toPlot <- toPlot[, colnames(ps$Y)[order(colSums(ps$Y), decreasing = TRUE)]]
    cov_order <- switch(as.character(whichPs), 
                        "3" = c("Genotype", "Population", "Plant.area", "Herbivory"), 
                        "4" = c("Genotype", "Population", "Plant.area", "Herbivory", "Random: Plant"), 
                        "5" = c("Population", "Plant.area", "Herbivory", "Random: Plant"),
                        "6" = c("Population", "Plant.area", "Herbivory", "Random: Plant", "Random: Genotype"),
                        "7" = c("Genotype", "Population", "Plant.area", "Herbivory", "Random: Plant"))
    toPlot <- toPlot[cov_order, ]
    round(rowMeans(toPlot*100), 2)

    wesandcols_sel <- wesandcols[1:nrow(toPlot)]
    pdf(file = file.path(dirs$fits, 
                         foldname, 
                         "figs",
                         paste0("varpart_", 
                                whichPs,
                                ".pdf")
                        ),
        bg = "white", 
        width = 7, 
        height = 5)
        par(family = "serif", mar = c(8,3,2,10), xpd = TRUE)
        barplot(toPlot, 
                col = wesandcols_sel,
                xpd = TRUE,
                las = 2)
        legend(x = 6.25, y = 1, 
               legend = rownames(toPlot)[nrow(toPlot):1],
               fill = wesandcols_sel[length(wesandcols_sel):1])
    dev.off()
}

# 3.3.2 residual correlations
library(circleplot)

if (whichPs == 7) {
    omgcors <- trap17:::computeAssociations_modified(ps)
    lev <- 1
    supportLevel <- 0
    toPlot <- ((omgcors[[lev]]$support > supportLevel)
              + (omgcors[[lev]]$support < (1-supportLevel)) > 0) * omgcors[[lev]]$mean
    toPlotDist <- list()
    for (i in 1:dim(toPlot)[3]) {
        toPlotDist[[i]] <- as.dist(toPlot[,,i])
        toPlotDist[[i]][which(toPlotDist[[i]] == 0)] <- NA
    }
} else {
    OmegaCor <- Hmsc:::computeAssociations(ps)
    lev <- 1
    supportLevel <- 0
    toPlot <- ((OmegaCor[[lev]]$support > supportLevel)
              + (OmegaCor[[lev]]$support < (1-supportLevel)) > 0) * OmegaCor[[lev]]$mean
    toPlotDist <- as.dist(toPlot)
    toPlotDist[which(toPlotDist == 0)] <- NA
}

if (whichPs == 7) {
    pdf(file = file.path(dirs$fits,
                         foldname,
                         "figs", 
                         paste0("omega_",
                         ps$rLNames[lev],
                         "_suppLev",
                         supportLevel * 100,
                         "_ps", 
                         whichPs,
                         ".pdf")),
        bg = "transparent", 
        width = 13, 
        height = 3)
        par(family = "serif", mfrow = c(1, length(toPlotDist)))
            for (i in 1:length(toPlotDist)) {
                circleplot(toPlotDist[prev_ord_genot][[i]], 
                           cluster = FALSE,
                           style = "classic",
                           plot.control = list(point.labels = TRUE,
                                               cex.point = 5,
                                               line.breaks = c(-1,0,1),
                                               line.cols = c("#00468b", "#c60032"),
                                               line.widths = 5))
            }
    dev.off()
} else {
    pdf(file = file.path(dirs$fits,
                         foldname,
                         "figs", 
                         paste0("omega_",
                         ps$rLNames[lev],
                         "_suppLev",
                         supportLevel * 100,
                         "_ps", 
                         whichPs,
                         ".pdf")),
        bg = "transparent", 
        width = 3, 
        height = 3)
        par(family = "serif")
        circleplot(toPlotDist, 
                   cluster = FALSE,
                   style = "classic",
                   plot.control = list(point.labels = TRUE,
                                       cex.point = 5,
                                       line.breaks = c(-1,0,1),
                                       line.cols = c("#00468b", "#c60032"),
                                       line.widths = 5))
    dev.off()
}

# 3.3.3 betas
summary(mpost$Beta)
betaPost_ps <- Hmsc:::getPostEstimate(ps, "Beta", q = c(0.05, 0.95))
betaMeans_ps <- betaPost_ps$mean
rownames(betaMeans_ps) <- colnames(ps$X)
supportLevel <- 0.75
betaSig <- ((betaPost_ps$support > supportLevel)
            + (betaPost_ps$support < (1-supportLevel)) > 0)
betaSig <- betaSig * 1            
rownames(betaSig) <- colnames(ps$X)

write.csv(round(betaMeans_ps, 2), 
          file = file.path(dirs$fits, 
                           foldname, 
                           "figs",
                            paste0("betameans_ps", whichPs, ".csv")))

betas <- round(betaMeans_ps * betaSig, 2)
betas    
betaPost_ps$q
colSums(dat$Y_pooled)

# 3.4 mixing
whichPs <- 3
ps <- pss[[whichPs]]

mpost <- Hmsc:::convertToCodaObject(ps)
par(family = "serif", mfrow = c(1, 1))
plot(mpost$Beta, auto.layout = FALSE, density = FALSE, ask = TRUE)

psfrs <- coda:::gelman.diag(mpost$Beta)$psrf
all((psfrs[,2] - psfrs[,1]) > 0)
sum(psfrs[,2] > 1.05)

#effectiveSize(mpost$Beta)
coda:::gelman.plot(mpost$Beta, ask = TRUE)

# 3.5 cooccurrence by genotype
library(cooccur)
cooc_test_genot <- list()
for (g in 1:length(unique(dat$X_pooled[,"Genotype"]))) {
    g1 <- sort(unique(dat$X_pooled[,"Genotype"]))[g]
    com <- t(dat$Y_pooled[which(dat$X_pooled[,"Genotype"]  == g1), ])
    tryCatch(
        expr = {
            cooc_test_genot[[g]] <- cooccur(mat = com,
                                            type = "spp_site",
                                            spp_names = TRUE,
                                            thresh = TRUE)
        },
        error = function(e) { print(e) },
        warning = function(w) { print(w) }
    )    
}
cooc_test <- cooccur(mat = t(dat$Y_pooled),
                     type = "spp_site",
                     spp_names = TRUE,
                     thresh = TRUE)


#cooc_test$positive
#cooc_test$results
#
cooc_test_genot[[2]]$results
cooc_test_genot[[2]]$positive
cooc_test_genot[[2]]$negative


# 3.6 raw data results figures
cooc_genot <- list()
for (g in 1:length(unique(dat$X_pooled[,"Genotype"]))) {
    g1 <- sort(unique(dat$X_pooled[,"Genotype"]))[g]
    cooc_genot[[g]] <- crossprod(dat$Y_pooled[which(dat$X_pooled[,"Genotype"]  == g1), ])
}

# NOTE: genotypes ordered by prevalence
pdf(file = file.path(dirs$wd,
                     "cooccurrences_genot.pdf"),
    bg = "transparent", 
    width = 13, 
    height = 3)
    par(family = "serif", mfrow = c(1, length(cooc_genot)))
    for (g in 1:length(cooc_genot)) {
        toPlotDist <- as.dist(cooc_genot[prev_ord_genot][[g]])
        circleplot:::circleplot(toPlotDist, 
                   cluster = FALSE,
                   style = "classic",
                   plot.control = list(point.labels = TRUE,
                                       cex.point = 5,
                                       line.breaks = c(-1,0,1,8,100),
                                       line.cols = c("#ffffff", 
                                                     "#faeaee", 
                                                     "#e795aa", 
                                                     "#c60032"),
                                       line.widths = 5))
    }
dev.off()

cooc_all <- crossprod(dat$Y_pooled)
toPlotDist <- as.dist(cooc_all)
pdf(file = file.path(dirs$wd,
                     "cooccurrences.pdf"),
    bg = "transparent", 
    width = 3, 
    height = 3)
    par(family = "serif")
        circleplot:::circleplot(toPlotDist, 
                   cluster = FALSE,
                   style = "classic",
                   plot.control = list(point.labels = TRUE,
                                       cex.point = 5,
                                       line.breaks = c(-1,0,1,8,100),
                                       line.cols = c("#ffffff", 
                                                     "#faeaee", 
                                                     "#e795aa", 
                                                     "#c60032"),
                                       line.widths = 5))
dev.off()

