# TRAP17 ANALYSIS 2020

# 1 DATA IMPORTATION
rm(list = ls(all = TRUE)) ; gc()
working_dir <- "/Users/anorberg/Documents/Zurich/UZH/TRAP/pkg/trap17-pkg"
setwd(working_dir)
library("trap17")
dirs <- set_dirs(working_dir = working_dir)
saveRDS(dirs, file = file.path(working_dir, "dirs.rds"))

# 1.1a process...
#dat <- process_data(filename = "TRAP17.csv",
#                    dirs = dirs,
#                    return_data = TRUE)

# 1.1b ...or load the processed data
?trapdata
dat <- trapdata
str(dat)

# 2 MODEL FITTING

# 2.1 sampling settings
sampling <- sampling_settings(totsamp = 1000,
                              trans = 500,
                              thn = 10,
                              nfolds = 5)
foldname <- create_name(study = "trap17",
                        totsamp = sampling$totsamp,
                        nfolds = sampling$nfolds, 
                        type = "fold")
create_directories(foldname = foldname, dirs = dirs)
saveRDS(sampling, 
        file = file.path(dirs$fits, foldname, "sampling.rds"))

# 2.2 model fitting and cross-validation
evals <- model_and_cv(dat = dat,
                      dirs = dirs,
                      variants = 1:5,
                      sampling = sampling,
                      returnCVs = TRUE,
                      saveCVs = TRUE) 
str(evals)
saveRDS(evals, 
        file = file.path(dirs$fits, foldname, "evals.rds"))

# 3 RESULTS
#rm(list = setdiff(ls(), c("working_dir", "dirs", "foldname", "sampling", "dat"))); gc()
#setwd(working_dir)
#library("trap17")
evals <- readRDS(file = file.path(file.path(dirs$fits, foldname), "evals.rds"))

# 3.1 cv-based R2s
tjurs <- lapply(lapply(evals, '[[', 1), '[[', 3)
cors <- lapply(lapply(evals, '[[', 2), colMeans, na.rm = TRUE)

# 3.2 Co-occurrence combinations
orig_combs <- co_occ_combs(partition = dat$X_pooled[,1:2],
                           Y_arr = array(dat$Y_pooled, 
                                         dim = c(nrow(dat$Y_pooled), 
                                                 ncol(dat$Y_pooled), 
                                                 1),
                                         dimnames = list(1:nrow(dat$Y_pooled),
                                                         colnames(dat$Y_pooled),
                                                         1))) 

modelled_combs <- vector(mode = "list", length = length(evals))

for (i in 1:length(evals)) {

    yarr <- evals[[i]]$cv_predictions
    dimnames(yarr) <- list(1:dim(yarr)[1],
                           colnames(dat$Y_pooled),
                           1:dim(yarr)[3])
    modelled_combs[[i]] <- co_occ_combs(partition = dat$X_pooled[,1:2],
                                        Y_arr = yarr) 
}

all_virus_combs <- as.character(unique(unlist(lapply(lapply(modelled_combs, 
                                                            dimnames), 
                                                            '[', 
                                                            1))))
all_virus_combs <- c("Empty", all_virus_combs[-which(all_virus_combs == "Empty")])

colrs <- load_colour_palette()
colmat <- cbind(all_virus_combs, colrs)
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
for (g in 1:4) {
    tmp11 <- c()
    for (p in 1:4) {
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
for (i in c(1:5)) {
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
    for (g in 1:4) {
        tmp11 <- c()
        for (p in 1:4) {
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
library(wesanderson)
wesandcols <- wes_palette("Cavalcanti1")[5:1]

pss <- load_objects_from_dir(path = dirs$fits, 
                             study = "trap17_totsamp1000",
                             obj_type = "ps")
names(pss)
whichPs <- 4
ps <- pss[[whichPs]]
    
# 3.3.1 variance partitioning
if (whichPs > 2) {
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
                col = wesandcols,
                xpd = TRUE,
                las = 2)
        legend(x = 6.25, y = 1, 
               legend = rownames(toPlot)[nrow(toPlot):1],
               fill = wesandcols[length(wesandcols):1])
    dev.off()
}

# 3.3.2 residual correlations
OmegaCor <- Hmsc:::computeAssociations(ps)
lev <- 1
OmegaCor[[lev]]$support
supportLevel <- 0
toPlot <- ((OmegaCor[[lev]]$support > supportLevel)
          + (OmegaCor[[lev]]$support < (1-supportLevel)) > 0) * OmegaCor[[lev]]$mean
toPlotDist <- as.dist(toPlot)
                        
pdf(file = file.path(dirs$fits,
                     foldname,
                     "figs", 
                     paste0("omega_suppLev",
                     supportLevel * 100,
                     "_ps", 
                     whichPs,
                     ".pdf")),
    bg = "transparent", 
    width = 4, 
    height = 4)
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

# 3.3.3 betas
betaPost_ps <- Hmsc:::getPostEstimate(ps, "Beta", q = c(0.1, 0.9))
betaMeans_ps <- betaPost_ps$mean
rownames(betaMeans_ps) <- colnames(ps$X)
supportLevel <- 0.7
betaSig <- ((betaPost_ps$support > supportLevel)
            + (betaPost_ps$support < (1-supportLevel)) > 0)
betaSig <- betaSig * 1            
rownames(betaSig) <- colnames(ps$X)
betas <- betaMeans_ps * betaSig
    

# 3.4 mixing
mpost <- Hmsc:::convertToCodaObject(ps)
names(mpost)
plot(mpost$Beta)
coda:::gelman.diag(mpost$Beta)$psrf
coda:::gelman.plot(mpost$Beta)
