# required libraries
library(AS.preprocessing); library(AS.angles); library(AS.circular); library(cluster)
setwd("~/Documents/ArchStats/Dissertation/sections/CS1-Genlis/img")

par(mar = c(2,2,0,0))

#=================================================================================================
# DATA CLEANING
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
save.features(final.classification, "Genlis-morph-final")
write.csv(centres, "Genlis-posthole-centres.csv", row.names = F)

#-------------------------------------------------------------------------------------------------
# reload features to avoid having to re-clean image later
genlis <- load.features("Genlis-morph-final")

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
vM.GoF(q.4, vm.mle$mu, vm.mle$kappa)
JP.GoF(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi)

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
# tests of similarity of quadrant distribution
# (no evidence that the two quadrants have different distributions)

q.samples <- list(q.4.a, q.4.b)
q.sizes <- c(length(q.4.a), length(q.4.b))

watson.common.mean.test(q.samples)                          # p = 0.81
wallraff.concentration.test(q.samples)                      # p = 0.74
mww.common.dist.LS(cs.unif.scores(q.samples), q.sizes)      # p = 0.88
watson.two.test(q.4.a, q.4.b)                               # p > 0.10
watson.two.test.rand(q.4.a, q.4.b, NR = 999)                # p = 0.49

JP.GoF(q.4.a, jp.mle$mu, jp.mle$kappa, jp.mle$psi)          # p > 0.15, p > 0.1
JP.GoF(q.4.b, jp.mle$mu, jp.mle$kappa, jp.mle$psi)          # p > 0.15, p > 0.1

#=================================================================================================
# GLOBAL VS LOCAL GRIDDING
#=================================================================================================
# fit uniform-von Mises mixture
em.u.vm <- EM.u.vonmises(q.4, k = 2)
plot.EM.vonmises(q.4, em.u.vm)

# proportion of distribution within certain range of pi
pmixt.vonmises <- function(q, model, p.from) {
    p <- rep(0, model$k)
    for (i in 1:length(model$k)) {
        p[i] <- model$alpha[i] * pvonmises(circular(q), circular(model$mu[i]), 
                                          model$kappa[i], from = circular(p.from), tol = 1e-06)
    }
    p
}
pmixt.vonmises(pi, em.u.vm, 0)

pmixedvonmises(pi, em.u.vm$mu[1], em.u.vm$mu[2], em.u.vm$kappa[1], em.u.vm$kappa[2], em.u.vm$alpha[1], from = 0)
# compare residuals
mvm.pp.res <- mvM.PP(q.4, em.u.vm$mu, em.u.vm$kappa, em.u.vm$alpha)

mean(mvm.pp.res^2); sd(mvm.pp.res^2)

#winner-takes-all clustering
em.clusts <- mvM.clusters(q.4, em.u.vm)

# use Euclidean distance to assess degree of spatial clustering
dist <- dist(pts)
em.sil <- silhouette(em.clusts, dist)
plot(em.sil, col = "black")

