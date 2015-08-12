
setwd("~/Documents/ArchStats/Dissertation/sections/gridding")
library(AS.circular); library(AS.angles)

# function to simulate points along the walls of a building, with some perturbation.
sim.building <- function(c1 = c(0,0), c2 = c(10,6), gap = 1, deg = 0, var = 0.1) {
    t <- max(c1[2], c2[2])
    b <- min(c1[2], c2[2])
    r <- max(c1[1], c2[1])
    l <- min(c1[1], c2[1])
    n.x <- (r-l)/gap
    n.y <- (t-b)/gap
    
    # get coordinates
    x <- c(l + c(0:n.x)*gap, rep(l, n.y-1),
           l + c(0:n.x)*gap, rep(r, n.y-1))
    y <- c(rep(b, n.x + 1), b + c(1:(n.y-1))*gap,
           rep(t, n.x + 1), b + c(1:(n.y-1))*gap)
    
    z <- cbind(x,y)
    z <- z + rnorm(length(z), 0, var)
    
    # apply rotation
    deg <- deg*pi/180
    rm <- matrix(c(cos(deg), sin(deg), -sin(deg), cos(deg)), nrow = 2, ncol = 2)
    z %*% rm
}

#=========================================================================
# EXAMPLE 1: ALL BUILDINGS ALIGNED, NO NOISE
set.seed(25247)
a <- sim.building(c(0,0), c(10,6), deg = 0)
b <- sim.building(c(-2,8), c(3,16), deg = 0)
c <- sim.building(c(11,4), c(15,15), deg = 0)
d <- sim.building(c(5,14), c(9,21), deg = 0)

z <- rbind(a,b,c,d)
#z <- rbind(z, cbind(runif(nrow(z)/2, floor(min(z[,1])), ceiling(max(z[,1]))),
#                    runif(nrow(z)/2, floor(min(z[,2])), ceiling(max(z[,2])))))

pdf(file = "sim-plot-1.pdf")
plot(z, pch = 20, asp = T, xlab = "", ylab = "")
dev.off()

k.1 <- k.nearest.angles(z, 1)
q <- circular(k.1[,-c(1:2)]) %% (2*pi)
pdf(file = "sim-q-plot-1.pdf")
plot(q, stack = T, sep = 0.05, pch = 20, shrink = 1.7)
lines(density.circular(q, bw = 15), lwd = 2, col = "red")
dev.off()

q.4 <- (4*q) %% (2*pi)
pdf(file = "sim-q4-plot-1.pdf")
plot(q.4, shrink = 1.7, stack = T, sep = 0.05, pch = 20)
lines(density.circular(q.4, bw = 15), lwd = 2, col = "red")
dev.off()

#======================================================================
# test for uniformity
rayleigh.test(q.4); kuiper.test(q.4); watson.test(q.4)

# test for symmetry
r.symm.test.stat(q.4); r.symm.test.boot(q.4, B = 999)

# obtain bias-corrected point estimates and confidence intervals
bc <- bc.ci.LS(q.4, alpha = 0.05)
(bc$mu[3] - bc$mu[1]) * 180 / pi            # range of mu in degrees
bc$beta2                                    # skewness
(bc$alpha2[1:3] - (bc$rho[c(1,3,2)]^4))     # kurtosis

# obtain ML estimates of von Mises parameters
# no bias correction, since sample size >= 16.
vm.mle <- mle.vonmises(q.4, bias = F, alp)
vm.mle$mu <- vm.mle$mu %% (2*pi)
q95 <- qnorm(0.975)
vm <- list(mu = c(est = vm.mle$mu, lower = (vm.mle$mu - q95*vm.mle$se.mu) %% (2*pi),
                  upper = (vm.mle$mu + q95*vm.mle$se.mu)) %% (2*pi),
           kappa = c(est = vm.mle$kappa, lower = vm.mle$kappa - q95*vm.mle$se.kappa, upper = vm.mle$kappa + q95*vm.mle$se.kappa))

jp.mle <- JP.mle(q.4)
jp <- JP.ci.nt(jp.mle, alpha = 0.05)

ests <- cbind(bc = rbind(bc$mu, bc$rho, A1inv(bc$rho), bc$beta2, bc$alpha2),
              vm = rbind(vm$mu, A1(vm$kappa), vm$kappa, NA, NA),
              jp = rbind(jp$mu %% (2*pi), A1(jp$kappa), jp$kappa, NA, jp$psi))
colnames(ests) <- c("est", "lower", "upper", "est.vm", "lower.vm", "upper.vm", "est.jp", "lower.jp", "upper.jp")
rownames(ests) <- c("$\\mu$", "$\\rho$", "$\\kappa$", "$\\beta_2$", "$\\alpha_2 / \\psi$\\footnote{$\\psi$ is given in place of $\\alpha_2$ for the Jones-Pewsey distribution}")
ests <- round(ests, 3)
write.csv(ests, file = "Simulated-ests.csv", row.names = T, quote = T)

# plots of density
circular.c.plot(q.4, l.pos = "bottomright")
Arrows(0, 0, 1.2 * 1.5 * cos(vm.mle$mu), 1.2 * 1.5 * 
           sin(vm.mle$mu), arr.type = "curved", lwd = 2, col = "blue")
Arrows(0, 0, 1.2 * 1.5 * cos(jp.mle$mu), 1.2 * 1.5 * 
           sin(jp.mle$mu), arr.type = "curved", lwd = 2, col = "red")

# centred linear plot
linear.c.plot((q.4 + (pi - bc$mu[1])) %% (2*pi))

vM.GoF(q.4, vm.mle$mu, vm.mle$kappa)
JP.GoF(q.4, jp.mle$mu, jp.mle$kappa, jp.mle$psi)
JP.GoF(q.4, jp.mle$mu, 0, 0)
vM.GoF(q.4, vm.mle$mu, 0)

AICc <- JP.psi.info(q.4, psi.0 = 0)$comparison[c(1,2,3,6),]
AICc[4,1] - AICc[4,2]

#======================================================================
# SPLIT DATA INTO QUADRANTS ACCORDING TO DIRECTION OF RAW ANGLE

mx <- circular(as.numeric(names(which.max(table(round(q, 1))))))
cutpoints <- sort(circular(mx + c(pi/4, 3*pi/4, 5*pi/4, 7*pi/4)) %% (2*pi))
quadrant <- rep(0, length(q))
quadrant[q > cutpoints[1] & q < cutpoints[2]] <- 1
quadrant[q > cutpoints[2] & q < cutpoints[3]] <- 2
quadrant[q > cutpoints[3] & q < cutpoints[4]] <- 3

q.4.a1 <- q.4[quadrant == 0]
q.4.b1 <- q.4[quadrant == 1]
q.4.a2 <- q.4[quadrant == 2]
q.4.b2 <- q.4[quadrant == 3]

q.4.a <- c(q.4.a1, q.4.a2)
q.4.b <- c(q.4.b1, q.4.b2)

pdf(file = "sim-quad-plot.pdf")
plot(q[quadrant == 0], pch = 20, stack = T, sep = 0.05, shrink = 2, bins = 90)
points(q[quadrant == 1], pch = 20, stack = T, sep = 0.05, shrink = 2, bins = 90, col = "blue")
points(q[quadrant == 2], pch = 20, stack = T, sep = 0.05, shrink = 2, bins = 90, col = "red")
points(q[quadrant == 3], pch = 20, stack = T, sep = 0.05, shrink = 2, bins = 90, col = "green4")
legend(1,-1, legend = c("Quadrant A", "Quadrant B", "Quadrant C", "Quadrant D"), col = c("Black", "Blue", "red", "green4"), pch = 20, bty = "n", cex = 1.4)
dev.off()

bc.a <- bc.ci.LS(q.4.a, alpha = 0.05)
bc.b <- bc.ci.LS(q.4.b, alpha = 0.05)

pdf(file = "sim-quad-plot-A.pdf")
circular.c.plot(q.4.a)
dev.off()

pdf(file = "sim-quad-plot-B.pdf")
circular.c.plot(q.4.b)
dev.off()

vm.mle.a <-  mle.vonmises(q.4.a, bias = T)
vm.mle.a$mu <- vm.mle.a$mu %% (2*pi)
jp.mle.a <- JP.mle(q.4.a)

bc.ci.LS(q.4.b, alpha = 0.05)
vm.mle.b <-  mle.vonmises(q.4.b, bias = T)
vm.mle.b$mu <- vm.mle.b$mu %% (2*pi)
jp.mle.b <- JP.mle(q.4.b)

#===============================================================
# tests of same distribution
watson.common.mean.test(list(q.4.a, q.4.b))
bc.a$mu[c(2,1,3)]
bc.b$mu[c(2,1,3)]

# some overlap between the two. Try estimating pooled mean as per Fisher.
m <- pooled.mean(list(q.4.a, q.4.b))
c(m$est - m$pm, m$est, m$est + m$pm) %% (2*pi)
bc$mu[c(2,1,3)]

# if bc$mu falls in the CI for the pooled mean, use bc$mu.
# Otherwise, use the pooled mean (global mean is being affected by other sectors)
