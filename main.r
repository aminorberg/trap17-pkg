
# 1 DATA: importation (and processing)
rm(list = ls(all = TRUE))
gc()

# define working directory
working_dir <- "/Users/anorberg/Documents/Zurich/UZH/TRAP/pkg/trap17-pkg"
setwd(working_dir)
library("trap17")

dirs <- set_dirs(working_dir = working_dir)
                 #raw_data = TRUE)
#saveRDS(dirs, file = file.path(working_dir, "dirs.rds"))

# 1.1a process...
#dat <- process_data(filename = "TRAP17.csv",
#                    dirs = dirs)

# 1.1b ...or load the processed data
dat <- trapdata
str(dat)
sum(rowSums(dat$Y_pooled) > 0) / nrow(dat$Y_pooled)
sum(rowSums(dat$Y_pooled) > 1) / sum(rowSums(dat$Y_pooled) > 0)

# 2 MODEL FITTING

# 2.1 sampling settings
# settings good for testing
sampling <- sampling_settings(totsamp = 100,
                              trans = 50,
                              thn = 1,
                              nchains = 1,
                              nfolds = 2)
# settings used for the submitted study
sampling <- sampling_settings(totsamp = 200000,
                             trans = 100000,
                             thn = 100,
                             nchains = 2,
                             nfolds = 10)

# 2.2 create a folder for the model with these settings and a figures folder under that
foldname <- create_name(study = "trap17",
                        totsamp = sampling$totsamp,
                        nfolds = sampling$nfolds, 
                        type = "fold")
create_directories(foldname = foldname, dirs = dirs)

# save also the sampling settings
saveRDS(sampling, 
        file = file.path(dirs$fits, foldname, "sampling.rds"))

# 2.3 model fitting and cross-validation
evals <- model_and_cv(dat = dat,
                      dirs = dirs,
                      variants = "ALL",
                      sampling = sampling,
                      saveCVs = TRUE)
saveRDS(evals, 
        file = file.path(dirs$fits, foldname, "evals.rds"))

# 3 RESULTS
foc_study <- "trap17_totsamp2e"
foldname <- "trap17_totsamp2e+05"

# 3.1.1 explanatory power
pss <- load_objects_from_dir(path = dirs$fits, 
                             study = foc_study,
                             obj_type = "ps")

eval_exp <- list()
for (i in 1:length(pss)) {
    preds_exp <- Hmsc::computePredictedValues(hM = pss[[i]], 
                                              alignPost = FALSE)
    eval_exp[[i]] <- Hmsc::evaluateModelFit(hM = pss[[i]], 
                                            predY = preds_exp)
}
eval_exp_arr <- lapply(eval_exp, simplify2array)
names(eval_exp_arr) <- paste0("ps", 1:ncol(eval_exp_means))
eval_exp_means <- simplify2array(lapply(eval_exp, lapply, mean))
colnames(eval_exp_means) <- paste0("ps", 1:ncol(eval_exp_means))
eval_exp_means

# 3.1.2 conditional predictions and predictive power
partition_sp <- 1:ncol(pss[[1]]$Y)
fit <- 2
vars <- trap17:::set_vars(study = "trap17",
                          fit = fit,
                          sampling = sampling)
cv_partition <- Hmsc:::createPartition(hM = pss[[2]], 
                                       nfolds = sampling$nfolds, 
                                       column = vars$partition)

if (!is.null(vars$covDepXvars)) {
    cv_preds <- trap17:::computePredictedValues_modified(hM = ps, 
                                                         partition = cv_partition,
                                                         partition.sp =  partition_sp,
                                                         expected = expected,
                                                         alignPost = FALSE)
} else {
    cv_preds <- Hmsc:::computePredictedValues(hM = pss[[fit]], 
                                              partition = cv_partition, 
                                              partition.sp = partition_sp,
                                              expected = expected,
                                              alignPost = TRUE)
}


# 3.1.3 cv-based R2s
evals <- readRDS(file = file.path(file.path(dirs$fits, foldname), "evals.rds"))
#evals <- evals[c(1:5,7)]
tjurs <- lapply(lapply(evals, '[[', 1), '[[', 3)
cors <- lapply(lapply(evals, '[[', 2), colMeans, na.rm = TRUE)
lapply(tjurs, mean)
lapply(cors, mean)


# 3.2 Co-occurrence combinations

# 3.2.1 Original combinations
orig_combs <- co_occ_combs(partition = dat$X_pooled[,c("Genotype","Population")],
                           Y_arr = array(dat$Y_pooled, 
                                         dim = c(nrow(dat$Y_pooled), 
                                                 ncol(dat$Y_pooled), 
                                                 1),
                                         dimnames = list(1:nrow(dat$Y_pooled),
                                                         colnames(dat$Y_pooled),
                                                         1))) 
saveRDS(orig_combs, 
        file = file.path(dirs$raw_data_figs, "orig_combs.rds"))

# 3.2.2 Predicted combinations
cv_preds <- load_objects_from_dir(path = dirs$fits, 
                                  study = foc_study,
                                  obj_type = "cv_preds")
names(cv_preds)
pss <- load_objects_from_dir(path = dirs$fits, 
                             study = foc_study,
                             obj_type = "ps")
names(pss)

# select the model you want to include in co-infection profile predictions
whichPs <- 3
cv_preds_variant <- cv_preds[whichPs]
str(cv_preds_variant)

# change names for the species for figures (*)
sp_new_nams <- c("Clo", "En", "Be", "Cap", "Cau") #(*)
sp_ord <- c("Clo", "Be", "Cap", "Cau", "En") #(*)

modelled_combs <- vector(mode = "list", length = length(cv_preds_variant))    
for (i in 1:length(modelled_combs)) {
    yarr <- cv_preds_variant[[i]]
    dimnames(yarr) <- list(1:dim(yarr)[1],
                           sp_new_nams,     #(*)
                           1:dim(yarr)[3])
    yarr <- yarr[,sp_ord,]                  #(*)
    modelled_combs[[i]] <- co_occ_combs(partition = dat$X_pooled[,c("Genotype",
                                                                    "Population")],
                                        Y_arr = yarr) 
}
saveRDS(modelled_combs, 
       file = file.path(dirs$fits, foldname, "modelled_combs.rds"))
modelled_combs <- readRDS(file = file.path(dirs$fits, foldname, "modelled_combs.rds"))

all_virus_combs <- as.character(unique(unlist(lapply(lapply(modelled_combs, 
                                                            dimnames), 
                                                            '[', 
                                                            1))))
all_virus_combs <- all_virus_combs[order(sapply(all_virus_combs, nchar), 
                                         decreasing = TRUE)]
all_virus_combs <- c("Empty", all_virus_combs[-which(all_virus_combs == "Empty")])

# 3.2.3 Plot original co-occurrence combinations (Fig 4)

# note: viruses, genotypes and populations are ordered by prevalence in the figures
sp_ord <- c("Clo", "Be", "Cap", "Cau", "En")
prev_ord_genot <- c(2, 4, 1, 3)
prev_ord_pop <- c(2, 3, 1, 4)
#

all_virus_combs <- c(all_virus_combs[-c(length(all_virus_combs):(length(all_virus_combs)-4))], 
                     sp_ord[length(sp_ord):1])

colrs <- load_colour_palette()
colmat <- cbind(all_virus_combs, colrs[[2]])
rownames(colmat) <- colmat[,1]
colmat_legend <- colmat

pdf(file = file.path(dirs$raw_data_figs, 
                     "orig_cooccs.pdf"),
    bg = "transparent", 
    width = 15, 
    height = 5)
    par(family = "serif", mfrow = c(1,4))
    for (g in prev_ord_genot) {
        tmp11 <- c()
        for (p in prev_ord_pop) {
            tmp12 <- orig_combs[, g, p][rownames(colmat)]
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

# 3.2.5 predicted co-occurrence combinations (Fig 4)
for (i in 1:length(modelled_combs)) {
    pdf(file = file.path(dirs$fits, 
                         foldname, 
                         "figs",
                         paste0("co_occs_ps", 
                                whichPs[i],
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

# 3.2.3 legend (Fig 4)
pdf(file = file.path(dirs$fits, 
                     foldname, 
                     "figs", 
                     "cooccurrencebars_legend.pdf"),
    bg = "transparent", 
    width = 5, 
    height = 10)
    par(family = "serif")
    plot(x = 1:10, 1:10, type = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("topleft", legend = colmat_legend[,1], fill = colmat_legend[,2], bty = 'n')
dev.off()

# 3.3 Results for the best model variant
pss <- load_objects_from_dir(path = dirs$fits, 
                             study = foc_study,
                             obj_type = "ps")
names(pss)
# select the model you want to include in co-infection profile predictions
whichPs <- 3
ps <- pss[[whichPs]]

library(wesanderson)
wesandcols <- wes_palette("Cavalcanti1")[5:1]

# 3.3.1 variance partitioning
if (whichPs > 2 & whichPs != 6) {
    sel <-  c(1, 2, rep(3, 3), rep(4, 3), 5) 
    selnames <- c("Intercept", 
                  "Plant.area", 
                  "Population", 
                  "Genotype", 
                  "Herbivory")
    if (whichPs == 5) {
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
                        "4" = c("Genotype", 
                                "Population", 
                                "Plant.area", 
                                "Herbivory", 
                                "Random: Plant"), 
                        "5" = c("Population", "Plant.area", "Herbivory", "Random: Plant"),
                        "6" = c("Genotype", 
                                "Population", 
                                "Plant.area", 
                                "Herbivory", 
                                "Random: Plant"))
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

if (whichPs == 2 | whichPs >= 4) {
    if (whichPs == 6) {
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

    if (whichPs == 6) {
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
} else {
    print(paste("Model variant", whichPs, "does not include latent variables"))
}

# 3.3.3 betas
betaPost_ps <- Hmsc:::getPostEstimate(ps, "Beta", q = c(0.05, 0.95))
betaMeans_ps <- betaPost_ps$mean
rownames(betaMeans_ps) <- colnames(ps$X)
supportLevel <- 0.95
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
write.csv(round(betas, 2), 
          file = file.path(dirs$fits, 
                           foldname, 
                           "figs",
                            paste0("betameans_", supportLevel, "_ps", whichPs, ".csv")))
   

# 3.3.4 mixing
whichPs <- 3
ps <- pss[[whichPs]]

mpost <- Hmsc:::convertToCodaObject(ps)
par(family = "serif", mfrow = c(1, 1))
plot(mpost$Beta, auto.layout = FALSE, density = FALSE, ask = TRUE)

psfrs <- coda:::gelman.diag(mpost$Beta)$psrf
# the lowest point estimate for the psfr  
min(psfrs[,1])
# highest 
max(psfrs[,1])
# all point estimates were below their corresponding upper confidence limits
all((psfrs[,2] - psfrs[,1]) > 0)

# 3.4 cooccurrence tests
library(cooccur)

# 3.4.1 cooccurrence test for full data
cooc_test <- cooccur(mat = t(dat$Y_pooled),
                     type = "spp_site",
                     spp_names = TRUE,
                     thresh = TRUE)
#cooc_test$results

# 3.4.2 cooccurrence by genotype test
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

# 3.4.3 cooccurrence by population test
cooc_test_pop <- list()
for (p in 1:length(unique(dat$X_pooled[,"Population"]))) {
    p1 <- sort(unique(dat$X_pooled[,"Population"]))[p]
    com <- t(dat$Y_pooled[which(dat$X_pooled[,"Population"]  == p1), ])
    tryCatch(
        expr = {
            cooc_test_pop[[p]] <- cooccur(mat = com,
                                            type = "spp_site",
                                            spp_names = TRUE,
                                            thresh = TRUE)
        },
        error = function(e) { print(e) },
        warning = function(w) { print(w) }
    )
}

# 3.5 raw data results figures

# 3.5.1 abundances by population and genotype (Fig 1)
#sp_ord <- c("Clo", "Be", "Cap", "Cau", "En")

abundances_genot <- list()
for (g in 1:length(unique(dat$X_pooled[,"Genotype"]))) {
    g1 <- sort(unique(dat$X_pooled[,"Genotype"]))[g]
    abundances_genot[[g]] <- colSums(dat$Y_pooled[which(dat$X_pooled[,"Genotype"]  == g1), ])
}
abundances_genot <- abundances_genot[prev_ord_genot]

abundances_pops <- list()
for (p in 1:length(unique(dat$X_pooled[,"Population"]))) {
    p1 <- sort(unique(dat$X_pooled[,"Population"]))[p]
    abundances_pops[[p]] <- colSums(dat$Y_pooled[which(dat$X_pooled[,"Population"]  == p1), ])
}
abundances_pops <- abundances_pops[prev_ord_pop]

colrs <- load_colour_palette()
colrs_rwplt <- colrs[[3]]
rownames(colrs_rwplt) <- colrs_rwplt[,1]

pdf(file = file.path(dirs$raw_data_figs, 
                     "abundances_by_pop.pdf"),
    bg = "transparent", 
    width = 4, 
    height = 5)
    par(family = "serif")
    barplot(simplify2array(abundances_pops),
            col = colrs_rwplt[,2],
            ylim = c(0, 200),
            xaxt = "n",
            yaxt = "n")
    axis(2, at = c(0, 50, 100, 200), labels = c("0", "50", "100", "320"), tick = TRUE, las = 2)
dev.off()

pdf(file = file.path(dirs$raw_data_figs, 
                     "abundances_by_genot.pdf"),
    bg = "transparent", 
    width = 4, 
    height = 5)
    par(family = "serif")
    barplot(simplify2array(abundances_genot),
            col = colrs_rwplt[,2],
            ylim = c(0, 200),
            xaxt = "n",
            yaxt = "n")
    axis(2, at = c(0, 50, 100, 200), labels = c("0", "50", "100", "320"), tick = TRUE, las = 2)
dev.off()

leg <- cbind(sp_ord, colrs_rwplt[,2])
leg <- leg[nrow(leg):1, ]
pdf(file = file.path(dirs$raw_data_figs, 
                     "abundances_legend.pdf"),
    bg = "transparent", 
    width = 3, 
    height = 5)
    par(family = "serif")
    plot(x = 1:10, 1:10, type = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("topleft", legend = c("Empty", leg[, 1]), fill = c("grey75", leg[, 2]), bty = 'n')
dev.off()


# 3.5.2 all co-occurrences (Fig 2)
library(circleplot)

pdf(file = file.path(dirs$raw_data_figs, 
                     "raw_coocc_legend.pdf"),
    bg = "transparent", 
    width = 5, 
    height = 3)
    par(family = "serif")
    plot(x = 1:10, 1:10, type = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("topleft", 
           legend = c("No co-occurrences",
                      "1 co-occurrence", 
                      "1-8 co-occurrences", 
                      ">8 co-occurrences"), 
           fill = c("#d9d9d9", 
                   "#faeaee", 
                   "#e795aa", 
                   "#c60032"), 
           bty = 'n')
dev.off()

cooc_all <- crossprod(dat$Y_pooled)
toPlotDist <- as.dist(cooc_all)

pdf(file = file.path(dirs$raw_data_figs,
                     "cooccurrences.pdf"),
    bg = "transparent", 
    width = 3, 
    height = 3)
    par(family = "serif")
        circleplot:::circleplot(toPlotDist, 
                   cluster = FALSE,
                   style = "classic",
                   plot.control = list(point.labels = TRUE,
                                       cex.point = 15,
                                       line.breaks = c(-1,0,1,8,100),
                                       line.cols = c("#d9d9d9", 
                                                     "#faeaee", 
                                                     "#e795aa", 
                                                     "#c60032"),
                                       line.widths = 5))
dev.off()

# 3.5.3 co-occurrences by genotype (Fig 2)

cooc_genot <- list()
for (g in 1:length(unique(dat$X_pooled[,"Genotype"]))) {
    g1 <- sort(unique(dat$X_pooled[,"Genotype"]))[g]
    cooc_genot[[g]] <- crossprod(dat$Y_pooled[which(dat$X_pooled[,"Genotype"]  == g1), ])
#    rownames(cooc_genot[[g]]) <- gsub("Ca", "Cau", rownames(cooc_genot[[g]]))
#    rownames(cooc_genot[[g]]) <- gsub("Pl", "Cap", rownames(cooc_genot[[g]]))
#    colnames(cooc_genot[[g]]) <- rownames(cooc_genot[[g]])
}
pdf(file = file.path(dirs$raw_data_figs,
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
                                       line.cols = c("#d9d9d9", 
                                                     "#faeaee", 
                                                     "#e795aa", 
                                                     "#c60032"),
                                       line.widths = 5))
    }
dev.off()

# 3.5.3 co-occurrences by population (Fig 2)

cooc_pop <- list()
for (p in 1:length(unique(dat$X_pooled[,"Population"]))) {
    p1 <- sort(unique(dat$X_pooled[,"Population"]))[p]
    cooc_pop[[p]] <- crossprod(dat$Y_pooled[which(dat$X_pooled[,"Population"]  == p1), ])
#    rownames(cooc_pop[[p]]) <- gsub("Ca", "Cau", rownames(cooc_pop[[p]]))
#    rownames(cooc_pop[[p]]) <- gsub("Pl", "Cap", rownames(cooc_pop[[p]]))
#    colnames(cooc_pop[[p]]) <- rownames(cooc_pop[[p]])
}
pdf(file = file.path(dirs$raw_data_figs,
                     "cooccurrences_pop.pdf"),
    bg = "transparent", 
    width = 13, 
    height = 3)
    par(family = "serif", mfrow = c(1, length(cooc_pop)))
    for (p in 1:length(cooc_pop)) {
        toPlotDist <- as.dist(cooc_pop[prev_ord_pop][[p]])
        circleplot:::circleplot(toPlotDist, 
                   cluster = FALSE,
                   style = "classic",
                   plot.control = list(point.labels = TRUE,
                                       cex.point = 5,
                                       line.breaks = c(-1,0,1,8,100),
                                       line.cols = c("#d9d9d9", 
                                                     "#faeaee", 
                                                     "#e795aa", 
                                                     "#c60032"),
                                       line.widths = 5))
    }
dev.off()



