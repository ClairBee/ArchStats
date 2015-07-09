
library(circular); library(FNN)

setwd("~/Documents/ArchStats/Mini-essays/Fitting a distribution")

org.par <- par()
set.cex.main <- 2
set.cex.lab <- 1.5
par(mar = c(0,0,0,0), cex.main = set.cex.main, cex.lab = set.cex.lab)

#============================================================================
# ANGULAR FUNCTIONS
# k.nearest.angles
# show.directions

# also import all functions in CIRCULAR FUNCTIONS

#============================================================================
# LOAD DATA & CONVERT MODULO PI/2
{
m <- read.csv("~/Documents/ArchStats/Fitting-dists/Genlis-15-07-01-details.csv")
p <- m[m$type == 0, 4:5]

k.1 <- k.nearest.angles(p, 1)

q <- circular(k.1[,-c(1,2)][!is.na(k.1[,-c(1,2)])])
q.4 <- (4*q) %% (2*pi)
n <- length(q.4)
}
#============================================================================
# get sample trigonometric moments; bias-corrected estimates; bootstrap estimates
{
bc.ests <- bc.point.estimates(q.4, sig = 0.05, symmetric = F)
vm.mle <-  mle.vonmises(q.4, bias = F)
vm.mle$mu <- vm.mle$mu %% (2*pi)
jp.mle <- JP.ci.nt(mle.jonespewsey(q.4), alpha = 0.05)

mu <- get.moments(q.4)$mu %% (pi*2)
r.bar <- get.moments(q.4)$r
beta.2 <- get.moments(q.4)$b2
alpha.2 <- get.moments(q.4)$a2

kappa <- vm.mle$kappa

vm.range <- (bc.ests$mu[3] - bc.ests$mu[2])*180/pi
jp.range <- (jp.mle$mu[3] - jp.mle$mu[2])*180/pi

# jp.mle.boot <- JP.ci.boot(q.4, alpha = 0.05, B = 9999)
# takes ages to run - values obtained are given below
    # mu:    3.313776  ( 3.258771,  3.368780)
    # kappa: 1.1871353 (0.8109857, 2.2702345)
    # psi:  -3.919586  (-5.393089, -2.862511)

jp.boot.range <- (jp.mle.boot$mu[3] - jp.mle.boot$mu[2])*180/pi

# show different CIs
plot(rep(0.7, 3), jp.mle.boot$kappa, pch = 20, ylim = c(-6,4), xlim = c(0.5, 1.5))
points(rep(0.758, 3), jp.mle$kappa, pch = 20, col = "red")
points(rep(1.2, 3), jp.mle.boot$psi, pch = 20)
points(rep(1.25, 3), jp.mle$psi, pch = 20, col = "red")
points(rep(0.975, 3), jp.mle.boot$mu, pch = 20)
points(rep(1.025, 3), jp.mle$mu, pch = 20, col = "red")
}
#============================================================================
# circular plot of angles modulo pi/2
{
pdf(file = "line-segment-plot.pdf")
show.directions(k.1)
dev.off()

# plot data, with densities
pdf(file = "plot-mod-pi-2.pdf")
plot(q.4, stack = T, sep = 0.05, pch = 20, shrink = 1.5, ylim = c(-1.2,0.8))
lines(density.circular(q.4, bw = 15), col = "blue")
arrows.circular(x = circular(bc.ests$mu[1]), shrink = 1.7)
arrows.circular(x = circular(jp.mle$mu[1]), shrink = 1.7, col = "cyan4")

curve.circular(dvonmises(x, mu = vm.mle$mu, kappa = vm.mle$kappa), n = 3600, add = T, lty = 2, col = "red", lwd = 2)
curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[1], psi = jp.mle$psi[1]), n = 3600, add = T, lty = 2, col = "cyan4", lwd = 2)
curve.circular(djonespewsey(x, mu = vm.mle$mu, kappa = vm.mle$kappa, psi = -1), n = 3600, add = T, lty = 1, col = "orange", lwd = 2)

legend("bottomright", cex = 0.8,
       legend = c("Kernel density estimate", "von Mises distribution", "Jones-Pewsey distribution", "Sample mean", expression(paste("Jones-Pewsey ", mu))),
       col = c("Blue", "Red", "cyan4", "Black", "cyan4"),
       lty = c(1, 2, 2, 1, 1),
       lwd = c(1, 2, 2, 1, 1))
dev.off()
}
#============================================================================
# reflective symmetry test
{
r.symm.test.stat(q.4)   # p = 0.6852
}

# uniformity tests
{
kuiper.test(q.4)
    # Test Statistic:  4.7332 
    # P-value < 0.01

watson.test(q.4)
    # Test Statistic: 1.8886 
    # P-value < 0.01

rao.spacing.test(q.4)
    # Test Statistic = 187.7545 
    # P-value < 0.001

rayleigh.test(q.4)
    # Test Statistic:  0.3392
    # P-value:  0
}
#============================================================================
# goodness-of-fit tests for von Mises distribution
{
watson.test(q.4, dist = "vonmises")
    # Test Statistic: 0.2818 
    # P-value < 0.01 

vm.unif <- circular(2*pi*pvonmises(q.4, circular(mu), kappa, from=circular(0), tol = 1e-06))
plot(vm.unif, stack = T, sep = 0.05, pch = 20, shrink = 1.5, ylim = c(-1.2,0.8))

kuiper.test(vm.unif)
    # Test Statistic:  2.0888 
    # P-value < 0.01

watson.test(vm.unif)
    # Test Statistic: 0.2823 
    # P-value < 0.01 

rao.spacing.test(vm.unif)
    # Test Statistic = 177.8316 
    # P-value < 0.001 

rayleigh.test(vm.unif)
    # Test Statistic: 0.0769 
    # P-value:  0.2017 
}
#============================================================================
# plot Jones-Pewsey with various values of psi
# based on max. likelihood estimation with various values of kappa & psi
{
plot(q.4, stack = T, sep = 0.05, pch = 20, shrink = 1.5, ylim = c(-1.3,0.7))
lines(density.circular(q.4, bw = 15), col = "blue")

curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[1], psi = jp.mle$psi[1]), add = T, lty = 1, col = "red", lwd = 2)
curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[1], psi = jp.mle$psi[3]), add = T, lty = 2, col = "darkred")
curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[1], psi = jp.mle$psi[2]), add = T, lty = 3, col = "darkred")

curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[2], psi = jp.mle$psi[3]), add = T, lty = 2, col = "orange")
curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[2], psi = jp.mle$psi[2]), add = T, lty = 3, col = "orange")

curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[3], psi = jp.mle$psi[3]), add = T, lty = 2, col = "green")
curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[3], psi = jp.mle$psi[2]), add = T, lty = 3, col = "green")

legend("bottomright", cex = 0.8,
       legend = c("Kernel density estimate", "MLE estimate",
                  expression(paste("MLE ", kappa, ", min ", psi)),
                  expression(paste("MLE ", kappa, ", max ", psi)),
                  expression(paste("Min ", kappa, ", min ", psi)),
                  expression(paste("Min ", kappa, ", max ", psi)),
                  expression(paste("Max ", kappa, ", min ", psi)),
                  expression(paste("Max ", kappa, ", max ", psi))),
       col = c("Blue", "Red", "darkred", "darkred", "orange", "orange", "green", "green"),
       lty = c(1, 1, 2, 3, 2, 3, 2, 3),
       lwd = c(2, 2, 1, 1, 1, 1, 1, 1))
}
# varying only psi to see effect on shape
{
plot(q.4, stack = T, sep = 0.05, pch = 20, shrink = 1.5, ylim = c(-1.3,0.7))
curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[1], psi = -1), add = T, lty = 1, col = "cornflowerblue", lwd = 2)
curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[1], psi = -2.5), add = T, lty = 1, col = "lightseagreen", lwd = 2)
curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[1], psi = -4), add = T, lty = 1, col = "chartreuse", lwd = 2)
curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[1], psi = -5.5), add = T, lty = 1, col = "goldenrod", lwd = 2)

legend("bottomright", cex = 0.8,
       legend = c(expression(paste(psi, " = -1")),
                  expression(paste(psi, " = -2.5")),
                  expression(paste(psi, " = -4")),
                  expression(paste(psi, " = -5.5"))),
       col = c("cornflowerblue", "lightseagreen", "chartreuse", "goldenrod"),
       lty = 1, lwd = 2)
}
#============================================================================
# P-P and Q-Q plots comparing vM & JP on same axis
{
edf <- ecdf(q.4)
ncon <- JP.NCon(jp.mle$kappa[1], jp.mle$psi[1]) 
jp.tdf <- 0 
jp.tqf <- 0
for (j in 1:length(q.4)) {
    jp.tdf[j] <- JP.df(q.4[j], jp.mle$mu[1], jp.mle$kappa[1], jp.mle$psi[1], ncon)
    jp.tqf[j] <- JP.qf(edf(q.4[j]), jp.mle$mu[1], jp.mle$kappa[1], jp.mle$psi[1], ncon)}

pdf(file = "P-P-plot.pdf")
plot(pvonmises(q.4, vm.mle$mu, vm.mle$kappa, from=circular(0), tol = 1e-06),
     edf(q.4), pch=20, xlim=c(0,1), ylim=c(0,1), asp = T, col = "grey",
     xlab = "Distribution function", ylab = "Empirical distribution function")
points(c(1:n)/(n+1), sort(q.4) / (2 * pi), pch = 20, col = "azure2")
points(jp.tdf, edf(q.4), pch=20)
lines(c(0,1), c(0,1), lwd=2, col = "red")
legend("bottomright", legend = c("Uniform", "von Mises", "Jones-Pewsey"), pch = 20, col = c("azure2", "grey", "black"))
dev.off()

pdf(file = "Q-Q-plot.pdf")
plot.default(qvonmises(edf(q.4), vm.mle$mu, vm.mle$kappa, from=circular(0), tol = 1e-06),
     q.4, pch=20, xlim=c(0,2*pi), ylim=c(0,2*pi), asp = T, col = "grey",
     xlab = "Quantile function", ylab = "Empirical quantile function")
points(jp.tqf, q.4, pch=20)
lines(c(0,2*pi), c(0,2*pi), lwd=2, col = "red")
legend("bottomright", legend = c("von Mises", "Jones-Pewsey"), pch = 20, col = c("grey", "black"))
dev.off()
}
#============================================================================
# goodness-of-fit tests for Jones-Pewsey distribution
{
JP.GoF.pvals(q.4, jp.mle$mu[1], jp.mle$kappa[1], jp.mle$psi[1])

# Kuiper            Test Statistic:  1.3658 
#                   P-value > 0.15 

# Watson            Test Statistic: 0.0746 
#                   P-value > 0.10 

# Rao spacing test  Test Statistic = 175.1859 
#                   P-value < 0.001  

# Rayleigh test     Test Statistic:  0.0428 
#                   P-value:  0.6085  
}
# could add in bootstrap test, but not really necessary
#============================================================================
# likelihood ratio tests & comparison of fit vs von Mises
{
JP.psi.LR.test(q.4, psi.0 = 0, alpha = 0.05)        # D = 28.929, p = 0
JP.psi.LR.test(q.4, psi.0 = -1, alpha = 0.05)       # D = 13.674, p = 0

JP.psi.info(q.4, psi.0 = 0)
JP.psi.info(q.4, psi.0 = -1)
    #                     AIC         BIC
    # Jones-Pewsey:   908.925     919.731
    # wrapped Cauchy: 920.599     927.803
    # von Mises:      935.854     943.058
}
#============================================================================
# try a Batschelet for comparison
Bat.mle <- Batmle(q.4)

# plot density
plot(q.4, stack = T, sep = 0.05, pch = 20, shrink = 1.5, ylim = c(-1.2,0.8))
lines(density.circular(q.4, bw = 15), col = "blue")
curve.circular(dvonmises(x, mu = vm.mle$mu, kappa = vm.mle$kappa), n = 3600, add = T, lty = 2, col = "chartreuse", lwd = 2)
curve.circular(djonespewsey(x, mu = circular(jp.mle$mu[1]), kappa = jp.mle$kappa[1], psi = jp.mle$psi[1]), n = 3600, add = T, lty = 2, col = "cyan4", lwd = 2)
curve.circular(BatPDF(x, Bat.mle$mu, Bat.mle$kappa, 1, ncon), add = T, n=3600, lty=2, lwd=2, col = "red")

bat.tdf <- 0; bat.tqf <- 0; bat.ncon <- BatNCon(Bat.mle$kappa, 1)
for (j in 1:length(q.4)) {
    bat.tdf[j] <- BatDF(q.4[j], Bat.mle$mu[1], Bat.mle$kappa[1], 1, bat.ncon)
    bat.tqf[j] <- BatQF(edf(q.4[j]), Bat.mle$mu, Bat.mle$kappa[1], 1, bat.ncon)}

# P-P & Q-Q plots
plot(pvonmises(q.4, vm.mle$mu, vm.mle$kappa, from=circular(0), tol = 1e-06),
     edf(q.4), pch=20, xlim=c(0,1), ylim=c(0,1), asp = T, col = "grey",
     xlab = "Distribution function", ylab = "Empirical distribution function")
points(bat.tdf, edf(q.4), pch=20, col = "cornflowerblue")
points(jp.tdf, edf(q.4), pch=20)
lines(c(0,1), c(0,1), lwd=2, col = "red")
legend("bottomright", legend = c("Uniform", "von Mises", "Jones-Pewsey", "Batschelet"), pch = 20, col = c("azure2", "grey", "black", "cornflowerblue"))

plot.default(qvonmises(edf(q.4), vm.mle$mu, vm.mle$kappa, from=circular(0), tol = 1e-06),
             q.4, pch=20, xlim=c(0,2*pi), ylim=c(0,2*pi), asp = T, col = "grey",
             xlab = "Quantile function", ylab = "Empirical quantile function")
points(bat.tqf, q.4, pch=20, col = "cornflowerblue")
points(jp.tqf, q.4, pch=20)
lines(c(0,2*pi), c(0,2*pi), lwd=2, col = "red")
legend("bottomright", legend = c("von Mises", "Jones-Pewsey", "Batschelet"), pch = 20, col = c("grey", "black", "cornflowerblue"))

# goodness-of-fit tests
BatGoF(q.4, Bat.mle$mu, Bat.mle$kappa, 1)
# Kuiper            Test Statistic:  0.9019 
#                   P-value > 0.15  

# Watson            Test Statistic: 0.0272 
#                   P-value > 0.10 

# Rao spacing test  Test Statistic = 173.419 
#                   P-value < 0.001 

# Rayleigh test     Test Statistic:  0.0053 
#                   P-value:  0.9924

BatCIBoot(q.4, B = 999)

# final parameter estimation (max. likelihood)
#                           mu              kappa                psi/nu
# von Mises         3.26 (3.05, 3.48)   0.72 (0.52, 0.93)           0
# Jones-Pewsey      3.31 (3.29, 3.34)   1.19 (0.59, 1.82)  -3.92 (-5.39, -2.45)
# Batschelet        3.26 (3.15, 3.37)   0.79 (0.64, 0.96)    1.70 (1.12, 2.19)