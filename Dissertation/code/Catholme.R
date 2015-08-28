# required libraries
library(AS.preprocessing); library(AS.angles); library(AS.circular)
setwd("~/Documents/ArchStats/Dissertation/img/CS2-Catholme")

par(mar = c(2,2,0,0))

#===========================================================================================
# DATA CLEANING
#===========================================================================================

# import map from JPEG image
catholme <- import.map("Catholme-cropped.jpg", threshold = 0.2, plot = T)

# exclude non-post-hole features
get.scale(catholme)                                   # identify scale marker
# no N-S marker on plan
catholme <- rescaled; remove(rescaled)

exclude.sparse.shapes(catholme, density = 0.55, lower = 3, plot = T)
get.postholes(sparse.shapes.classified)

# no further cleaning necessary

catholme <- final.classification
save.features(catholme, "Catholme-final")
catholme <- load.features("Catholme-final")
get.postholes(catholme)
#---------------------------------------------------------------------------------------

# further cleaning: exclude isolated and non-feature points
dist.filter <- filter.by.distance(centres)
nn.filter <- filter.by.2nn(centres)


# extract angles
k.1 <- k.nearest.angles(centres[dist.filter & nn.filter,], 1)
q <- circular(k.1[,-c(1,2)][!is.na(k.1[,-c(1,2)])]) %% (2*pi)
q.4 <- (4*q) %% (2*pi)              # convert axial to circular data by 'wrapping'

#====================================================================================
# TESTS TO FIT AND SELECT A GLOBAL MODEL
#====================================================================================

# test for uniformity and symmetry
rayleigh.test(q.4); kuiper.test(q.4); watson.test(q.4)

r.symm.test.stat(q.4); r.symm.test.boot(q.4, B = 999)

#------------------------------------------------------------------------------------
# parameter estimation

bc <- bc.ci.LS(q.4, alpha = 0.05)

(bc$mu[3] - bc$mu[1]) * 180 / pi            # range of mu in degrees
bc$beta2                                    # skewness
(bc$alpha2[1:3] - (bc$rho[c(1,3,2)]^4))

vm.mle <- mle.vonmises(q.4, bias = F)
vm.mle$mu <- vm.mle$mu %% (2*pi)

# calculate confidence interval using normal theory
q95 <- qnorm(0.975)
vm <- list(mu = c(est = vm.mle$mu,
                  lower = (vm.mle$mu - q95*vm.mle$se.mu) %% (2*pi),
                  upper = (vm.mle$mu + q95*vm.mle$se.mu)) %% (2*pi),
           kappa = c(est = vm.mle$kappa,
                     lower = vm.mle$kappa - q95*vm.mle$se.kappa,
                     upper = vm.mle$kappa + q95*vm.mle$se.kappa))

jp.mle <- JP.mle(q.4)
jp <- JP.ci.nt(jp.mle, alpha = 0.05)

#------------------------------------------------------------------------------------
# Goodness-of-fit tests
vM.GoF(q.4, vm.mle$mu, vm.mle$kappa)
JP.GoF(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi)

AICc <- JP.psi.info(q.4, psi.0 = 0)$comparison[c(1,2,3,6),]
AICc[4,1] - AICc[4,2]

