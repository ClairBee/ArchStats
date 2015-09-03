# required libraries
library(AS.preprocessing); library(AS.angles); library(AS.circular)
setwd("~/Documents/ArchStats/Dissertation/sections/CS2-Catholme/img")

par(mar = c(2,2,0,0))

#=================================================================================================
# DATA CLEANING
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

get.postholes(sparse.shapes.classified)
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

# get normalising constant for max.likelihood distribution
jp.ncon <- JP.NCon(jp.mle$kappa, jp.mle$psi)

#-------------------------------------------------------------------------------------------------
# Goodness-of-fit tests
vM.GoF(q.4, vm.mle$mu, vm.mle$kappa)
JP.GoF(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi)



# plot and calculate residuals
vm.pp.res <- vM.PP(q.4.adj, vm.mle$mu, vm.mle$kappa)
jp.pp.res <- JP.PP(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi)
vm.qq.res <- vM.QQ(q.4, vm.mle$mu, vm.mle$kappa)
jp.qq.res <- JP.QQ(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi)

# mean squared error and standard deviation
mean(vm.pp.res^2); sd(vm.pp.res^2)
mean(jp.pp.res^2); sd(jp.pp.res^2)
mean(vm.qq.res^2); sd(vm.qq.res^2)
mean(jp.qq.res^2); sd(jp.qq.res^2)



