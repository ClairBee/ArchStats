# required libraries
library(AS.preprocessing); library(AS.angles); library(AS.circular); library(fpc)
setwd("~/Documents/ArchStats/Dissertation/sections/CS2-Catholme/img")

par(mar = c(2,2,0,0))

#=================================================================================================
# DATA CLEANING FOR CATHOLME PLAN
#=================================================================================================

# import map from JPEG image
catholme <- import.map("Catholme-cropped.jpg", threshold = 0.2, plot = F)

# these functions do not overwrite the original site object
# exclude non-post-hole features
get.scale(catholme)                                 # identify scale marker
# no N-S marker on plan
catholme <- rescaled; remove(rescaled)

exclude.sparse.shapes(catholme, density = 0.55, lower = 3, plot = F)
# no further cleaning necessary
catholme <- sparse.shapes.classified

get.postholes(catholme)
save.features(final.classification, "Catholme-final")

#-------------------------------------------------------------------------------------------------
# reload features to avoid having to re-clean image later
catholme <- load.features("Catholme-final")
get.postholes(catholme)

# further cleaning: exclude isolated points
dist.filter <- filter.by.distance(centres)
pts <- centres[dist.filter,]

# extract angles
k.1 <- k.nearest.angles(pts, 1)
q <- circular(k.1[,-c(1,2)][!is.na(k.1[,-c(1,2)])]) %% (2*pi)
q.4 <- (4*q) %% (2*pi)              # convert axial to circular data by 'wrapping'


#=================================================================================================
# TESTS TO FIT AND SELECT MODELS
#=================================================================================================
# test for uniformity and symmetry
rayleigh.test(q.4)                  # p = 0
kuiper.test(q.4)                    # p < 0.01
watson.test(q.4)                    # p < 0.01

r.symm.test.stat(q.4)               # p = 0.019

#-------------------------------------------------------------------------------------------------
# parameter estimation
bc <- bc.ci.LS(q.4, alpha = 0.05)

(bc$mu[3] - bc$mu[1]) * 180 / pi            # range of mu in degrees: +- 12
bc$beta2                                    # skewness
(bc$alpha2[1:3] - (bc$rho[c(1,3,2)]^4))     # non-zero excess kurtosis

vm.mle <- mle.vonmises(q.4, bias = F)
vm.mle$mu <- vm.mle$mu %% (2*pi)
q95 <- qnorm(0.975)
vm <- list(mu = c(est = vm.mle$mu,
                  lower = (vm.mle$mu - q95*vm.mle$se.mu) %% (2*pi),
                  upper = (vm.mle$mu + q95*vm.mle$se.mu)) %% (2*pi),
           kappa = c(est = vm.mle$kappa,
                     lower = vm.mle$kappa - q95*vm.mle$se.kappa,
                     upper = vm.mle$kappa + q95*vm.mle$se.kappa))

jp.mle <- JP.mle(q.4)
jp <- JP.ci.nt(jp.mle, alpha = 0.05)

# normalising constant for max.likelihood distribution
jp.ncon <- JP.NCon(jp.mle$kappa, jp.mle$psi)

#-------------------------------------------------------------------------------------------------
# Goodness-of-fit tests
vM.GoF.boot(q.4, vm.mle$mu, vm.mle$kappa, B = 999)                          # p = 0.001, 0.001
JP.GoF.boot(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi, B = 999)              # p = 0.013, 0.040

#=================================================================================================
# DIVIDE SITE INTO SMALLER REGIONS AND RE-TEST: ANY EVIDENCE OF GRIDDING?
#=================================================================================================
# Density-based clustering
db.clust <- dbscan(pts, MinPts = 4, eps = 5)$cluster
db.clusters <- sort(table(db.clust[db.clust > 0]), decreasing = T)
db.clusters <- data.frame("id" = names(db.clusters)[db.clusters >= 25],
                          "size" = db.clusters[db.clusters >= 25])

#-------------------------------------------------------------------------------------------------
# create table of cluster results
for (i in 1:nrow(db.clusters)) {
    id <- rownames(db.clusters)[i]
    a <- q.4[db.clust == id]

    vm.a <- mle.vonmises(a, bias = F)
    jp.a <- JP.mle(a)
    bc.a <- bc.sample.statistics(a)
        
    # tests of uniformity and reflective symmetry
    if (rayleigh.test(a)$p.val < 0.05) {db.clusters$unif.rayl[i] <- "Non-U"} else
                                       {db.clusters$unif.rayl[i] <- "Uniform"}
    if (kuiper.test(a)$statistic > 1.747) {db.clusters$unif.kuip[i] <- "Non-U"} else
                                          {db.clusters$unif.kuip[i] <- "Uniform"}
    if (watson.test(a)$statistic > 0.187) {db.clusters$unif.wats[i] <- "Non-U"} else
                                          {db.clusters$unif.wats[i] <- "Uniform"}
    
    if (r.symm.test.stat(a)$p.val < 0.05) {db.clusters$symm[i] <- "Skewed"} else
                                          {db.clusters$symm[i] <- "Symmetric"}
        
    if ("Uniform" %in% db.clusters[i,3:5] | "Skewed" %in% db.clusters[i,6]) {
        db.clusters$vM.fit.k[i] <-NA;  db.clusters$vM.fit.w[i] <- NA
        db.clusters$JP.fit.k[i] <- NA; db.clusters$JP.fit.w[i] <- NA
        db.clusters$vM[i] <- NA; db.clusters$JP[i] <- NA
        db.clusters$mu[i] <- NA; db.clusters$rho[i] <- NA
    } else {
        db.clusters$mu[i] <- bc.a$mu %% (2*pi)
        db.clusters$rho[i] <- bc.a$rho
        vm.fit <- vM.GoF.boot(a, vm.a$mu, vm.a$kappa, B = 999, show.progress = T)
        jp.fit <- JP.GoF.boot(a, jp.a$mu, jp.a$kappa, jp.a$psi, B = 999, show.progress = T)
        
        db.clusters$vM.fit.k[i] <- vm.fit[1]; db.clusters$vM.fit.w[i] <- vm.fit[2]
        if (max(db.clusters[i,9:10]) < 0.05) {
            db.clusters$vM[i] <- NA
        } else {
            db.clusters$vM[i] <- paste("(", round(vm.a$mu %% (2*pi), 2), ", ",
                                       round(vm.a$kappa, 2), ")", sep = "")
        }
        
        db.clusters$JP.fit.k[i] <- jp.fit[1]; db.clusters$JP.fit.w[i] <- jp.fit[2]
        if (max(db.clusters[i,12:13]) < 0.05) {
            db.clusters$JP[i] <- NA
        } else {
            db.clusters$JP[i] <- paste("(", round(jp.a$mu %% (2*pi), 2), ", ",
                                       round(jp.a$kappa, 2), ", ", round(jp.a$psi, 2), 
                                       ")", sep = "") 
        }
    }
}
write.csv(db.clusters, "DBclust-tests.csv", row.names = F, quote = T)

q.samples <- list(); q.sizes <- c()
for (i in 1:length(cand)) {
    q.samples[[i]] <- q.4[db.clust == cand[i]]
    q.sizes[i] <- length(q.samples[[i]])
}

#=================================================================================================
# TEST PERPENDICULARITY
#=================================================================================================
# split points into quadrants
mx <- circular(as.numeric(names(which.max(table(round(q, 1))))))

cutpoints <- sort(circular(mx + c(pi/4, 3*pi/4, 5*pi/4, 7*pi/4)) %% (2*pi))
quadrant <- rep(0, length(q))
quadrant[findInterval(q, cutpoints) %in% c(1,3)] <- 1

# test perpendicularity
quad.tests <- matrix(ncol = 8, nrow = length(q.samples),
                     dimnames = list(q.sizes, c("dir", "conc", "ray.a", "wats.a", "kuip.a", "ray.b", "wats.b", "kuip.b")))
for (i in 1:length(q.samples)) {
    qa <- q.4[quadrant == 0 & db.clust == cand[i]]
    qb <- q.4[quadrant == 1 & db.clust == cand[i]]
    quad.tests[i,1] <- signif(watson.common.mean.test(list(qa, qb))$p.val, 2)
    quad.tests[i,2] <- signif(wallraff.concentration.test(list(qa, qb))$p.val, 2)
    
    if (rayleigh.test(qa)$p.val < 0.05) {quad.tests[i,3] <- "-"} else {quad.tests[i,3] <- "Uniform"}
    if (watson.test(qa)$statistic > 0.187) {quad.tests[i,4] <- "-"} else {quad.tests[i,4] <- "Uniform"}
    if (kuiper.test(qa)$statistic > 1.747) {quad.tests[i,5] <- "-"} else {quad.tests[i,5] <- "Uniform"}
    
    if (rayleigh.test(qb)$p.val < 0.05) {quad.tests[i,6] <- "-"} else {quad.tests[i,6] <- "Uniform"}
    if (watson.test(qb)$statistic > 0.187) {quad.tests[i,7] <- "-"} else {quad.tests[i,7] <- "Uniform"}
    if (kuiper.test(qb)$statistic > 1.747) {quad.tests[i,8] <- "-"} else {quad.tests[i,8] <- "Uniform"}
}

quad.tests
