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
JP.GoF.boot(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi, B = 999)              # p = 0.022, 0.046

#=================================================================================================
# DIVIDE SITE INTO SMALLER REGIONS AND RE-TEST: ANY EVIDENCE OF GRIDDING?
#=================================================================================================
# Density-based clustering
db.clust <- dbscan(pts, MinPts = 4, eps = 5)$cluster
db.clusters <- sort(table(db.clust[db.clust > 0]), decreasing = T)
db.clusters <- data.frame("id" = names(db.clusters)[db.clusters >= 25],
                          "size" = db.clusters[db.clusters >= 25])

#-------------------------------------------------------------------------------------------------
# create table of cluster results (including per-quadrant)
for (i in 1:nrow(db.clusters)) {
    id <- rownames(db.clusters)[i]
    a <- q.4[db.clust == id]
    
    # tests of uniformity and reflective symmetry
    db.clusters$unif.rayl[i] <- rayleigh.test(a)$p.val
    if (kuiper.test(a)$statistic > 1.747) {db.clusters$unif.kuip[i] <- "Non-U"} else
                                          {db.clusters$unif.kuip[i] <- "Uniform"}
    if (watson.test(a)$statistic > 0.187) {db.clusters$unif.wats[i] <- "Non-U"} else
                                          {db.clusters$unif.wats[i] <- "Uniform"}
    
    if (length(a) < 50) {db.clusters$symm[i] <- r.symm.test.boot(a)$p.val} else
                        {db.clusters$symm[i] <- r.symm.test.stat(a)$p.val}
        
    # identify quadrants for this particular data set
    mx <- circular(as.numeric(names(which.max(table(round(q[db.clust == id], 1))))))
    cutpoints <- sort(circular(mx + c(pi/4, 3*pi/4, 5*pi/4, 7*pi/4)) %% (2*pi))

    qa <- a[!(findInterval(q[db.clust == id], cutpoints) %in% c(1,3))]
    qb <- a[findInterval(q[db.clust == id], cutpoints) %in% c(1,3)]
    
    db.clusters$rayl.a[i] <- rayleigh.test(qa)$p.val
    if (kuiper.test(qa)$statistic > 1.747) {db.clusters$kuip.a[i] <- "Non-U"} else
                                           {db.clusters$kuip.a[i] <- "Uniform"}
    if (watson.test(qa)$statistic > 1.747) {db.clusters$wats.a[i] <- "Non-U"} else
                                           {db.clusters$wats.a[i] <- "Uniform"}
    
    db.clusters$rayl.b[i] <- rayleigh.test(qb)$p.val
    if (kuiper.test(qb)$statistic > 1.747) {db.clusters$kuip.b[i] <- "Non-U"} else
                                           {db.clusters$kuip.b[i] <- "Uniform"}
    if (watson.test(qb)$statistic > 0.187) {db.clusters$wats.b[i] <- "Non-U"} else
                                           {db.clusters$wats.b[i] <- "Uniform"}
    db.clusters$mean[i] <- signif(watson.common.mean.test(list(qa, qb))$p.val, 2)
    db.clusters$conc[i] <- signif(wallraff.concentration.test(list(qa, qb))$p.val, 2)
}
write.csv(db.clusters, "DBclust-tests.csv", row.names = F, quote = T)
