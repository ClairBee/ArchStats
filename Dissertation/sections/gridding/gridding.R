
setwd("~/Documents/ArchStats/Dissertation/sections/gridding")

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

verts <- rbind(c(1,2), c(6,11), c(11,4), c(15,15))

set.seed(24747)
a <- sim.building(c(0,0), c(10,6), deg = 0)
b <- sim.building(c(-2,8), c(3,16), deg = 0)
c <- sim.building(c(11,4), c(15,15), deg = 0)
d <- sim.building(c(5,14), c(9,21), deg = 5)

z <- rbind(a,b,c,d)
z <- rbind(z, cbind(runif(nrow(z)/2, floor(min(z[,1])), ceiling(max(z[,1]))),
                    runif(nrow(z)/2, floor(min(z[,2])), ceiling(max(z[,2])))))

pdf(file = "sim-plot.pdf")
plot(z, pch = 20, asp = T, xlab = "", ylab = "")
dev.off()

k.1 <- k.nearest.angles(z, 1)
q <- circular(k.1[,-c(1:2)])
pdf(file = "sim-q-plot.pdf")
plot(q, stack = T, sep = 0.05, pch = 20, shrink = 1.7)
lines(density.circular(q, bw = 15), lwd = 2, col = "red")
dev.off()

q.4 <- (4*q) %% (2*pi)
pdf(file = "sim-q4-plot.pdf")
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
bc$rho^4

(bc$alpha2 - (bc$rho^4)) / ((1-bc$rho)^2)
