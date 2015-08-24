# required libraries
library(AS.preprocessing); library(AS.angles); library(AS.circular)
setwd("~/Documents/ArchStats/Dissertation/img")

par(mar = c(2,2,0,0))

#====================================================================================
# DATA CLEANING
#====================================================================================

# ADD IN ALL PARAMETERS ----------------------------------------------------------------

# import map from JPEG image
genlis <- import.map("Genlis-cropped.jpg", threshold = 0.2, plot = F)

# exclude non-post-hole features
get.scale(genlis)                                   # identify scale marker
genlis.NS <- get.NS.axis(rescaled)                  # identify N-S marker
exclude.sparse.shapes(NS.marked)
remove.annotations(sparse.shapes.classified)
remove.tall.features(annotations.removed)
extend.annotations(tall.features.removed)

# convert to points
get.postholes(with.extensions)
genlis <- final.classification

# further cleaning: exclude isolated and non-feature points
dist.filter <- filter.by.distance(centres)
nn.filter <- filter.by.2nn(centres)

# extract angles
# NEED TO UPDATE FUNCTIONS TO GIVE CORRECT OUTPUT FORMAT -------------------------------
k.1 <- k.nearest.angles(centres[dist.filter $ nn.filter,], 1)
q <- circular(k.1[,-c(1,2)][!is.na(k.1[,-c(1,2)])]) %% (2*pi)
q.4 <- (4*q) %% (2*pi)              # convert axial to circular data by 'wrapping'


#====================================================================================
# TESTS TO FIT AND SELECT MODELS
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

vm.mle <- mle.vonmises(q.4, bias = F, alp)
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

#------------------------------------------------------------------------------------
# export results to .csv
ests <- cbind(bc = rbind(bc$mu, bc$rho, A1inv(bc$rho), bc$beta2, bc$alpha2),
              vm = rbind(vm$mu, A1(vm$kappa), vm$kappa, NA, NA),
              jp = rbind(jp$mu %% (2*pi), A1(jp$kappa), jp$kappa, NA, jp$psi))
colnames(ests) <- c("est", "lower", "upper", 
                    "est.vm", "lower.vm", "upper.vm",
                    "est.jp", "lower.jp", "upper.jp")
rownames(ests) <- c("mu", "rho", "kappa", "beta2", "alpha2", "psi")
ests <- round(ests, 3)
write.csv(ests, file = "Genlis-simulated-ests.csv", row.names = T, quote = T)

#------------------------------------------------------------------------------------
# Goodness-of-fit tests
vM.GoF(q.4, vm.mle$mu, vm.mle$kappa)
JP.GoF(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi)

AICc <- JP.psi.info(q.4, psi.0 = 0)$comparison[c(1,2,3,6),]
AICc[4,1] - AICc[4,2]


#====================================================================================
# LINEARITY VS PERPENDICULARITY
#====================================================================================

# split points into quadrants
# find approximate modal angle
mx <- circular(as.numeric(names(which.max(table(round(q, 1))))))

cutpoints <- circular(mx + c(pi/4, 3*pi/4, 5*pi/4, 7*pi/4)) %% (2*pi)
quadrant <- rep(0, length(q))
quadrant[q > cutpoints[1] & q < cutpoints[2]] <- 1
quadrant[q > cutpoints[3] & q < cutpoints[4]] <- 1

q.4.a <- q.4[quadrant == 0]
q.4.b <- q.4[quadrant == 1]

#-----------------------------------------------------------------------------------
# get bias-corrected and ML point estimates for each quarter
bc.a <- bc.sample.statistics(q.4.a, symmetric = F)
vm.a <-  mle.vonmises(q.4.a, bias = T)
vm.a$mu <- vm.a$mu %% (2*pi)
jp.a <- JP.mle(q.4.a)

bc.b <- bc.sample.statistics(q.4.b, symmetric = F)
vm.b <-  mle.vonmises(q.4.b, bias = T)
vm.b$mu <- vm.b$mu %% (2*pi)
jp.b <- JP.mle(q.4.b)

#-----------------------------------------------------------------------------------
# tests of similarity of quadrant distribution
# (no evidence that the two quadrants have different distributions)

q.samples <- list(q.4.a, q.4.b)
q.sizes <- c(length(q.4.a), length(q.4.b))

watson.common.mean.test(q.samples)                          # p = 0.76
wallraff.concentration.test(q.samples)                      # p = 0.77
mww.common.dist.LS(cs.unif.scores(q.samples), q.sizes)      # p = 0.58
watson.two.test(q.4.a, q.4.b)                               # p > 0.10
watson.two.test.rand(q.4.a, q.4.b, NR = 999)                # p = 0.68


#====================================================================================
# GLOBAL VS LOCAL GRIDDING
#====================================================================================
# divide data into regions using midpoint method

# use pca-based method to cluster into metres
# (could then try this with points from linear features added in)

# cluster points into grid of 1-2m across (param based on density?)
# and test pairwise similarity of each grid section to get a clustering?
# plot clusters as in pca analysis.

# change starting point and grid size slightly?
