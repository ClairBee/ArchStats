# required libraries
library(AS.preprocessing); library(AS.angles); library(AS.circular)
library(fpc)
setwd("~/Documents/ArchStats/Dissertation/sections/CS1-Genlis/img")

par(mar = c(2,2,0,0))

#=================================================================================================
# DATA CLEANING FOR GENLIS PLAN
#=================================================================================================

# import map from JPEG image
genlis <- import.map("Genlis-cropped.jpg", threshold = 0.2, plot = F)

# these functions do not overwrite the original site object
# exclude non-post-hole features
get.scale(genlis)                                   # identify scale marker
genlis.NS <- get.NS.axis(rescaled)                  # identify N-S marker
genlis <- NS.marked; remove(rescaled, NS.marked)

exclude.sparse.shapes(genlis, density = 0.55, lower = 3, plot = F)
genlis.closing <- feature.closing(sparse.shapes.classified, plot.progress = F)
fill.broken.boundary(features.closed, s = 0.2, plot.progress = F)
genlis <- boundary.filled

get.postholes(genlis)
save.features(genlis, "Genlis-final")
write.csv(centres, "Genlis-posthole-centres.csv", row.names = F)

#-------------------------------------------------------------------------------------------------
# reload features to avoid having to re-clean image later
genlis <- load.features("Genlis-final")

# preserves row names, rather than loading centres from csv
get.postholes(genlis)

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

r.symm.test.stat(q.4)               # p = 0.374

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

# get normalising constant for max.likelihood distribution
jp.ncon <- JP.NCon(jp.mle$kappa, jp.mle$psi)

#-------------------------------------------------------------------------------------------------
# Goodness-of-fit tests
vM.GoF.boot(q.4, vm.mle$mu, vm.mle$kappa, B = 9999)                          # p = 0.001, 0.001
JP.GoF.boot(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi, B = 9999)              # p = 0.259, 0.161 

# plot and calculate residuals
vm.pp.res <- vM.PP(q.4, vm.mle$mu, vm.mle$kappa)
jp.pp.res <- JP.PP(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi)
vm.qq.res <- vM.QQ(q.4, vm.mle$mu, vm.mle$kappa)
jp.qq.res <- JP.QQ(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi)

# mean squared error and standard deviation
mean(vm.pp.res^2); sd(vm.pp.res^2)
mean(jp.pp.res^2); sd(jp.pp.res^2)
mean(vm.qq.res^2); sd(vm.qq.res^2)
mean(jp.qq.res^2); sd(jp.qq.res^2)


#=================================================================================================
# LINEARITY VS PERPENDICULARITY
#=================================================================================================
# split points into quadrants
# find approximate modal angle
mx <- circular(as.numeric(names(which.max(table(round(q, 1))))))

cutpoints <- sort(circular(mx + c(pi/4, 3*pi/4, 5*pi/4, 7*pi/4)) %% (2*pi))
quadrant <- rep(0, length(q))
quadrant[findInterval(q, cutpoints) %in% c(1,3)] <- 1

q.4.a <- q.4[quadrant == 0]
q.4.b <- q.4[quadrant == 1]

#-------------------------------------------------------------------------------------------------
# get bias-corrected and ML point estimates for each quarter
bc.a <- bc.sample.statistics(q.4.a)
vm.a <-  mle.vonmises(q.4.a, bias = T)
vm.a$mu <- vm.a$mu %% (2*pi)
jp.a <- JP.mle(q.4.a)
jp.ci.a <- JP.ci.nt(jp.a, alpha = 0.05)

bc.b <- bc.sample.statistics(q.4.b)
vm.b <-  mle.vonmises(q.4.b, bias = T)
vm.b$mu <- vm.b$mu %% (2*pi)
jp.b <- JP.mle(q.4.b)
jp.ci.b <- JP.ci.nt(jp.b, alpha = 0.05)
#-------------------------------------------------------------------------------------------------
# tests of uniformity
rayleigh.test(q.4.a); rayleigh.test(q.4.b)              # p = 0         p = 0
kuiper.test(q.4.a); kuiper.test(q.4.b)                  # p < 0.01      p < 0.01
watson.test(q.4.a); watson.test(q.4.b)                  # p < 0.01      p < 0.01

r.symm.test.stat(q.4.a); r.symm.test.stat(q.4.b)        # p = 0.542     p = 0.475


# tests of similarity of quadrant distribution
# (no evidence that the two quadrants have different distributions)
q.samples <- list(q.4.a, q.4.b); q.sizes <- c(length(q.4.a), length(q.4.b))

watson.common.mean.test(q.samples)                          # p = 0.81
wallraff.concentration.test(q.samples)                      # p = 0.74
mww.common.dist.LS(cs.unif.scores(q.samples), q.sizes)      # p = 0.88
watson.two.test(q.4.a, q.4.b)                               # p > 0.10
watson.two.test.rand(q.4.a, q.4.b, NR = 999)                # p = 0.49

JP.GoF.boot(q.4.a, jp.mle$mu, jp.mle$kappa, jp.mle$psi, B = 999)          # p = 0.122, 0.077
JP.GoF.boot(q.4.b, jp.mle$mu, jp.mle$kappa, jp.mle$psi, B = 999)          # p = 0.229, 0.188 

#=================================================================================================
# FIT MIXTURE MODEL
#=================================================================================================
# fit uniform-von Mises mixture
em.vm <- EM.vonmises(q.4, k = 2)
em.u.vm <- EM.u.vonmises(q.4, k = 2)
plot.EM.vonmises(q.4, em.u.vm)

# compare residuals
uvm.pp.res <- mvM.PP(q.4, em.u.vm$mu, em.u.vm$kappa, em.u.vm$alpha)
mvm.pp.res <- mvM.PP(q.4, em.vm$mu, em.vm$kappa, em.vm$alpha)

mean(uvm.pp.res^2); sd(uvm.pp.res^2)
mean(mvm.pp.res^2); sd(mvm.pp.res^2)
mean(jp.pp.res^2); sd(jp.pp.res^2)

# AIC for each model
n <- length(q.4)
(2*3) - (2 * jp.mle$maxll) + (2*3*4)/(n-3-1)       # k = 3; AICc = 808.4157
(2*5) - (2 * em.vm$log.lh) + (2*5*6)/(n-5-1)       # k = 5; AICc = 804.0046
(2*3) - (2 * em.u.vm$log.lh) + (2*3*4)/(n-3-1)     # k = 3; AICc = 800.8829

# winner-takes-all clustering based on uniform-von Mises mixture
em.clusts <- mvM.clusters(q.4, em.u.vm)

#-------------------------------------------------------------------------------------------------
# test von Mises component for perpendicularity
clust.a <- q.4[quadrant == 0 & em.clusts == 2]
clust.b <- q.4[quadrant == 1 & em.clusts == 2]

c.samples <- list(clust.a, clust.b); c.sizes <- c(length(clust.a), length(clust.b))

watson.common.mean.test(c.samples)                          # p = 0.494
wallraff.concentration.test(c.samples)                      # p = 0.044
mww.common.dist.LS(cs.unif.scores(c.samples), c.sizes)      # p = 0.072
watson.two.test(clust.a, clust.b)                           # 0.05 < p < 0.10
watson.two.test.rand(clust.a, clust.b, NR = 999)            # p = 0.053

# check confidence intervals for rho
bc.vm.a <- bc.ci.LS(clust.a, alpha = 0.05)
bc.vm.b <- bc.ci.LS(clust.b, alpha = 0.05)
bc.vm.a; bc.vm.b

#=================================================================================================
# ASSESS DEGREE OF GRIDDING
#=================================================================================================
# plot clustered points on grid
plot(pts[em.clusts == 1,], pch = 20, col = "grey")
points(pts[em.clusts == 2,], pch = 20, col = "black")

# spatial clustering using DBscan
db.clust <- dbscan(pts, MinPts = 4, eps = 4.65)$cluster
points(pts, col = db.clust + 1)

xt <- xtabs(~., data = cbind("component" = c("uniform", "von Mises")[em.clusts],
                             "region" = db.clust))
sweep(xt, 2, colSums(xt), "/")

#             region
# component           0         1         2         3         4         5
# uniform     0.8571429 0.4545455 0.5000000 0.4021739 0.5000000 0.5000000
# von Mises   0.1428571 0.5454545 0.5000000 0.5978261 0.5000000 0.5000000

#-------------------------------------------------------------------------------------------------
# subset into density-based clusters and test for common orientation
q.db1 <- q.4[db.clust == 1]
q.db3 <- q.4[db.clust == 3]

db.samples <- list(q.db1, q.db3); db.sizes <- c(length(q.db1), length(q.db3))

watson.common.mean.test(db.samples)                          # p = 0.213
wallraff.concentration.test(db.samples)                      # p = 0.176
mww.common.dist.LS(cs.unif.scores(db.samples), db.sizes)     # p = 0.546
watson.two.test(q.db1, q.db3)                                # p > 0.10
watson.two.test.rand(q.db1, q.db3, NR = 999)                 # p = 0.497
