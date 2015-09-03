setwd("~/Documents/ArchStats/Dissertation/sections/CS2-Catholme/img")
# Row numbers refer to Catholme.R
point.colour <- "grey"; JP.colour = "red"; vM.colour = "blue"; bins <- 90; BW = 30
#------------------------------------------------------------------------------------
# 33: plot centres
{
pdfheight <- round(ymax(catholme$features) / 10, 0)/2
pdfwidth <- round(xmax(catholme$features) / 10, 0)/2

pdf("Catholme-postholes.pdf", height = pdfheight, width = pdfwidth)
plot(catholme$features, col = "cornflowerblue", cex.axis = 1.3, asp = F, legend = F, frame.plot = F)
points(centres, pch = 20)
points(centres[!dist.filter,], col = "red2", pch = 20, cex = 1.1, lwd = 2)
legend("topright", legend = c("Post-hole", "Excluded by distance filter"), bty = "n",
       pch = 20, col = c("black", "red2"), cex = 1.3)
dev.off()
}

#------------------------------------------------------------------------------------
# 73: circular and linear plots of the transformed data, with MLE distributions
# centre data at mu
q.4.adj <- q.4 + (as.numeric(q.4 < (bc$mu[1] - pi)) * 2*pi)
{
# bin data
cuts <- c(0:bins) * 2 * pi/bins
b <- circular(((cuts[1:bins] + cuts[2:(bins + 1)])/2)[findInterval(q.4, cuts)])
b.l <- matrix(b)
kd <- cbind(density.circular(q.4, bw = BW)$x,
            density.circular(q.4, bw = BW)$y)

pdf(file = "Q-circ-plot.pdf")
plot(circular(((cuts[1:bins] + cuts[2:(bins + 1)])/2)[findInterval(q, cuts)]),
     axes = F, shrink = 2.5, col = point.colour, stack = T, sep = 0.05, ylim = c(-1.2,0.8))
axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.2,
              labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                         expression(paste("3", pi, "/2"))))
Arrows(0.7 * cos(mx + c(0, pi/2, pi, 3*pi/2)), 0.7 * sin(mx + c(0, pi/2, pi, 3*pi/2)),
       0.8 * cos(mx + c(0, pi/2, pi, 3*pi/2)), 0.8 * sin(mx + c(0, pi/2, pi, 3*pi/2)), lty = 2, col = "seagreen")
lines(density.circular(q, bw = BW), lwd = 3)
legend("bottom", bty = "n", cex = 1.3, col = c("black", "seagreen"), lty = c(1,1,4), lwd = c(3, 2),
       legend = c("Kernel density estimate", expression(paste("Modal direction + 0, ", pi, "/2, ", pi, ", 3", pi, "/2"))))
dev.off()

pdf(file = "Q4-circ-plot.pdf")
plot(q.4, axes = F, shrink = 2.5, col = point.colour, stack = T, sep = 0.05, ylim = c(-1.2,0.8))
axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.2,
              labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                         expression(paste("3", pi, "/2"))))
lines(density.circular(q.4, bw = BW), lwd = 3)
curve.circular(dvonmises(x, mu = vm.mle$mu, kappa = vm.mle$kappa), n = 3600, add = T,
               lty = 2, col = vM.colour, lwd = 3)
curve.circular(djonespewsey(x, mu = jp.mle$mu, kappa = jp.mle$kappa, psi = jp.mle$psi),
               n = 3600, add = T, lty = 4, col = JP.colour, lwd = 3)
legend(-1.5,-2, bty = "n", cex = 1.3, col = c("black", vM.colour, JP.colour), lty = c(1,2,4), lwd = 3,
       legend = c("Kernel density estimate", "von Mises distribution", "Jones-Pewsey distribution"))
dev.off()

# adjust data to centre the histogram on the sample mean direction
q.4.adj <- b.l + ((b.l < (bc$mu[1] - pi)) * 2*pi)
kd2 <- cbind(x = c(kd[,1], kd[,1] + (2*pi)),
             y = rep(kd[,2], 2))
kd3 <- kd2[kd2[,1] > (min(q.4.adj) - pi/40) & kd2[,1] < (max(q.4.adj) + pi/40),]

xlabl <- c(expression(paste("-", pi, "/2")), 0, 
           expression(paste(pi, "/2")), expression(paste(pi)),
           expression(paste("3", pi, "/2")), expression(paste("2", pi)),
           expression(paste("5", pi, "/2")), expression(paste("3", pi)))
xbreaks <- c(-0.5, 0, 0.5, 1,1.5,2,2.5, 3) * pi

pdf(file = "Q4-linear-plot.pdf")
hist(matrix(q.4.adj), xaxt = "none", ylim = c(0,0.5), col = point.colour, breaks = 40, cex.axis = 1.5,
     xlab = "Transformed angle (radians)", border = "darkgrey", cex.lab = 1.5, main = "", freq = F)
lines(kd3, col = "black", lwd = 3)
curve(dvonmises(x, mu = vm.mle$mu, kappa = vm.mle$kappa), n = 3600, add = T, lty = 2, 
      col = vM.colour, lwd = 3)
curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, 
      lty = 4, col = JP.colour, lwd = 3)
axis(1, at = xbreaks, cex.axis = 1.5, labels = xlabl)
#legend("topright", bty = "n", cex = 1.3, col = c("black", vM.colour, JP.colour), lty = c(1,2,4), lwd = 3,
#       legend = c("Kernel density estimate", "von Mises distribution", "Jones-Pewsey distribution"))
dev.off()
}

#------------------------------------------------------------------------------------
# 86: probability plots
{
    n <- length(q.4)
    edf <- ecdf(q.4)
    ncon <- JP.NCon(jp.mle$kappa, jp.mle$psi)
    
    vm.tdf <- pvonmises(q.4, vm.mle$mu, vm.mle$kappa, from = circular(0), tol = 1e-06)
    jp.tdf <- 0
    for (j in 1:n) {
        jp.tdf[j] <- JP.df(q.4[j], jp.mle$mu, jp.mle$kappa, jp.mle$psi, ncon)
    }
    
    vm.tqf <- qvonmises(edf(q.4), vm.mle$mu, vm.mle$kappa, from = circular(0), tol = 1e-06)
    jp.tqf <- 0
    for (j in 1:n) {
        jp.tqf[j] <- JP.qf(edf(q.4)[j], jp.mle$mu, jp.mle$kappa, jp.mle$psi, ncon)
    }
    
    hist(matrix(q.4), xaxt = "none", ylim = c(0,0.5), col = point.colour, breaks = 40, cex.axis = 1.5,
         xlab = "Transformed angle (radians)", border = "darkgrey", cex.lab = 1.5, main = "", freq = F)
    lines(kd2, col = "black", lwd = 3)
    curve(dvonmises(x, mu = vm.mle$mu, kappa = vm.mle$kappa), n = 3600, add = T, lty = 2, 
          col = vM.colour, lwd = 3)
    curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, 
          lty = 4, col = JP.colour, lwd = 3)
    axis(1, at = xbreaks, cex.axis = 1.5, labels = xlabl)

}




#------------------------------------------------------------------------------------
cl <- cbind(catholme$feature.types[,1], NA)
cl[cl[,1] %in% rownames(pts),2] <- q.4
r.clusters <- reclassify(catholme$features, cl)

moran.global <- Moran(r.clusters)

moran <- c()
for (i in 1:100) {
    cl[cl[,1] %in% rownames(pts),2] <- rvonmises(nrow(pts), vm.mle$mu, vm.mle$kappa)
    moran[i] <- Moran(reclassify(catholme$features, cl))
}


