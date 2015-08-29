setwd("~/Documents/ArchStats/Dissertation/sections/CS2-Catholme/img")
# Row numbers refer to Genlis.R
point.colour <- "grey"; JP.colour = "red"; vM.colour = "blue"; bins <- 90; BW = 15
#------------------------------------------------------------------------------------
# 33: plot centres
pdfheight <- round(ymax(catholme$features) / 10, 0)/2
pdfwidth <- round(xmax(catholme$features) / 10, 0)/2

pdf("Catholme-postholes.pdf", height = pdfheight, width = pdfwidth)
par(mar = c(2,2,0,0))
plot(catholme$features, col = "white", cex.axis = 1.3, asp = F, legend = F, frame.plot = F)
points(centres[dist.filter & nn.filter,], pch = 20)
dev.off()

#------------------------------------------------------------------------------------
# 73: export results to .csv
ests <- cbind(bc = rbind(bc$mu, bc$rho, A1inv(bc$rho), bc$beta2, bc$alpha2, NA),
              vm = rbind(vm$mu, A1(vm$kappa), vm$kappa, NA, NA, NA),
              jp = rbind(jp$mu %% (2*pi), A1(jp$kappa), jp$kappa, NA, NA, jp$psi))
colnames(ests) <- c("est", "lower", "upper", 
                    "est.vm", "lower.vm", "upper.vm",
                    "est.jp", "lower.jp", "upper.jp")
rownames(ests) <- c("mu", "rho", "kappa", "beta2", "alpha2", "psi")
ests <- round(ests, 3)
write.csv(ests, file = "../../csv/Catholme-ests.csv", row.names = T, quote = T)

#------------------------------------------------------------------------------------
# 73: circular and linear plots of the transformed data, with MLE distributions

# bin data
cuts <- c(0:bins) * 2 * pi/bins
b <- circular(((cuts[1:bins] + cuts[2:(bins + 1)])/2)[findInterval(q.4, cuts)])
b.l <- matrix(b)
kd <- cbind(density.circular(q.4, bw = BW)$x,
            density.circular(q.4, bw = BW)$y)
kd2 <- cbind(x = c(kd[,1], kd[,1] + (2*pi)),
             y = rep(kd[,2], 2))


pdf(file = "Q4-circ-plot.pdf")
plot(b, axes = F, shrink = 2, col = point.colour, stack = T, sep = 0.05, ylim = c(-1.4,0.6))
axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.2,
              labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                         expression(paste("3", pi, "/2"))))
lines(density.circular(q.4, bw = BW), lwd = 3)
curve.circular(dvonmises(x, mu = vm.mle$mu, kappa = vm.mle$kappa), n = 3600, add = T,
               lty = 2, col = vM.colour, lwd = 3)
curve.circular(djonespewsey(x, mu = jp.mle$mu, kappa = jp.mle$kappa, psi = jp.mle$psi),
               n = 3600, add = T, lty = 4, col = JP.colour, lwd = 3)
legend("bottom", bty = "n", cex = 1.3, col = c("black", vM.colour, JP.colour), lty = c(1,2,4), lwd = 3,
       legend = c("Kernel density estimate", "von Mises distribution", "Jones-Pewsey distribution"))
dev.off()

# adjust data to centre the histogram on the sample mean direction
q.4.adj <- b.l + ((b.l < (bc$mu[1] - pi)) * 2*pi)
kd3 <- kd2[kd2[,1] > (min(q.4.adj) - pi/40) & kd2[,1] < (max(q.4.adj) + pi/40),]

pdf(file = "Q4-linear-plot.pdf")
hist(matrix(q.4.adj), xaxt = "none", col = point.colour, breaks = 40, cex.axis = 1.5,
     xlab = "Transformed angle (radians)", border = "darkgrey", cex.lab = 1.5, main = "", freq = F)
lines(kd3, col = "black", lwd = 3)
curve(dvonmises(x, mu = vm.mle$mu, kappa = vm.mle$kappa), n = 3600, add = T, lty = 2, 
      col = vM.colour, lwd = 3)
curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, 
      lty = 4, col = JP.colour, lwd = 3)
axis(1, at = c(0,0.5,1,1.5,2,2.5,3,3.5,4) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi, " = 0")),
                expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
#legend("topright", bty = "n", cex = 1.3, col = c("black", vM.colour, JP.colour), lty = c(1,2,4), lwd = 3,
#       legend = c("Kernel density estimate", "von Mises distribution", "Jones-Pewsey distribution"))
dev.off()

pdf(file = "Q-circ-plot.pdf")
plot(circular(((cuts[1:bins] + cuts[2:(bins + 1)])/2)[findInterval(q, cuts)]),
     axes = F, shrink = 2, col = point.colour, stack = T, sep = 0.05, ylim = c(-1.4,0.6))
axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.2,
              labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                         expression(paste("3", pi, "/2"))))
lines(density.circular(q, bw = BW), lwd = 3)
legend("bottom", bty = "n", cex = 1.3, col = c("black"), lty = c(1,2,4), lwd = 3,
       legend = c("Kernel density estimate"))
dev.off()