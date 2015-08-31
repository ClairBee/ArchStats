setwd("~/Documents/ArchStats/Dissertation/sections/CS1-Genlis/img")
# Row numbers refer to Genlis.R
point.colour <- "grey"; JP.colour = "red"; vM.colour = "blue"; bins <- 90; BW = 15

#--------------------------------------------------------------------------------------------
pdfheight <- round(ymax(genlis$features) / 10, 0)
pdfwidth <- round(xmax(genlis$features) / 10, 0)

plot.to.pdf <- function(site, pdf.filename, h, w, show.points = F) {
    site$feature.types[site$feature.types[,2] %in% c(2,3), 2] <- 4
    cols <- c("black", "black", "red", "red", "cornflowerblue", "cornflowerblue", "red")
    desc <- c("Unclassified", "Post-hole", "Scale / N-S axis", "OBSOLETE", "Removed", "Removed previously", "Changed")
    ind <- sort(unique(site$feature.types[, 2])) + 1
    l.pch <- 20
    
    pdf(file = pdf.filename, height = h, width = w)
    plot(reclassify(site$features, rcl = site$feature.types), col = cols[ind], cex.axis = 1.3, asp = F, legend = F, frame.plot = F)
    if (show.points) {
        xy <- data.frame(cbind(id = getValues(site$features), xyFromCell(site$features, 1:ncell(site$features))))
        mids <- ddply(xy[!is.na(xy$id), ], .(id), summarise, xm = mean(x), ym = mean(y))
        points(mids[site$feature.types[,2] == 4,2:3], cex = 1.2, lwd = 2,
               col = cols[7], xlim = c(0, xmax(site$features)), ylim = c(0, ymax(site$features)))
        ind[ind == 6] <- 7
        l.pch <- c(rep(20, length(ind)-1),1)
    }
    legend("topleft", legend = desc[ind], col = cols[ind], ncol = 2, pch = l.pch, cex = 1.3, bty = "n", text.font = 3)
    dev.off()
}

#--------------------------------------------------------------------------------------------
# 20: plot sparse shapes classified
plot.to.pdf(sparse.shapes.classified, "Genlis-sparse.pdf", h = pdfheight, w = pdfwidth)

# 23: plot features identified as closed features
ft <- features.closed$feature.types
ft[sparse.shapes.classified$feature.types[,2] %in% c(2,3,4),2] <- 5
plot.to.pdf(list(features = genlis$features, feature.types = ft), "Genlis-after-closing.pdf", h = pdfheight, w = pdfwidth, show.points = T)

# 24: plot with broken boundary filled
ft <- boundary.filled$feature.types
ft[features.closed$feature.types[,2] %in% c(2,3,4),2] <- 5
plot.to.pdf(list(features = genlis$features, feature.types = ft), "Genlis-boundary-filled.pdf", h = pdfheight, w = pdfwidth, show.points = T)

#28: post-hole centres
pdf("Genlis-1-postholes.pdf", height = pdfheight, width = pdfwidth)
plot(genlis$features, col = "cornflowerblue", cex.axis = 1.3, asp = F, legend = F, frame.plot = F)
points(centres, pch = 20)
#points(centres[!nn.filter,], col  = "purple", lwd = 2, pch = 4)
points(centres[!dist.filter,], col = "red", lwd = 2)
# points(centres[!nn.filter,], pch = 4, col = "red", lwd = 2)
legend("topleft", legend = c("Post-hole", "Excluded by distance filter", "Exluded by angular filter"), bty = "n",
       pch = c(20, 1), col = c("black", "red"), cex = 1.3)
dev.off()

#------------------------------------------------------------------------------------
# 87: export results to .csv
ests <- cbind(bc = rbind(bc$mu, bc$rho, A1inv(bc$rho), bc$beta2, bc$alpha2, NA),
              vm = rbind(vm$mu, A1(vm$kappa), vm$kappa, NA, NA, NA),
              jp = rbind(jp$mu %% (2*pi), A1(jp$kappa), jp$kappa, NA, NA, jp$psi))
colnames(ests) <- c("est", "lower", "upper", 
                    "est.vm", "lower.vm", "upper.vm",
                    "est.jp", "lower.jp", "upper.jp")
rownames(ests) <- c("mu", "rho", "kappa", "beta2", "alpha2", "psi")
ests <- round(ests, 3)
write.csv(ests, file = "../../csv/Genlis-ests.csv", row.names = T, quote = T)
#------------------------------------------------------------------------------------

# 73: circular and linear plots of the transformed data, with MLE distributions
{
# bin data
cuts <- c(0:bins) * 2 * pi/bins
b <- circular(((cuts[1:bins] + cuts[2:(bins + 1)])/2)[findInterval(q.4, cuts)])
kd <- cbind(density.circular(q.4, bw = BW)$x,
            density.circular(q.4, bw = BW)$y)
kd2 <- cbind(x = c(kd[,1], kd[,1] + (2*pi)),
             y = rep(kd[,2], 2))


pdf(file = "Q4-circ-plot.pdf")
plot(b, axes = F, shrink = 2, col = point.colour, stack = T, sep = 0.05, ylim = c(-1.2,0.8))
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

# no need to centre the histogram on the sample mean direction
pdf(file = "Q4-linear-plot.pdf")
hist(matrix(b), xaxt = "none", col = point.colour, breaks = 40, cex.axis = 1.5,
     xlab = "Transformed angle (radians)", border = "darkgrey", cex.lab = 1.5, main = "", freq = F)
lines(kd2, col = "black", lwd = 3)
curve(dvonmises(x, mu = vm.mle$mu, kappa = vm.mle$kappa), n = 3600, add = T, lty = 2, 
      col = vM.colour, lwd = 3)
curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, 
      lty = 4, col = JP.colour, lwd = 3)
axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
dev.off()

pdf(file = "Q-circ-plot.pdf")
plot(circular(((cuts[1:bins] + cuts[2:(bins + 1)])/2)[findInterval(q, cuts)]),
     axes = F, shrink = 2, col = point.colour, stack = T, sep = 0.05, ylim = c(-1.1,0.9))
axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.2,
              labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                         expression(paste("3", pi, "/2"))))
lines(density.circular(q, bw = BW), lwd = 3)
legend("bottom", bty = "n", cex = 1.3, col = c("black"), lty = c(1,2,4), lwd = 2,
       legend = c("Kernel density estimate"))
dev.off()
}


#===========================================================================================
#===========================================================================================
# plots of von Mises & Jones-Pewsey fit
n <- length(q.4)
edf <- ecdf(data)

vm.tdf <- pvonmises(q.4, vm.mle$mu, vm.mle$kappa, from = circular(0), tol = 1e-06)
jp.tdf <- 0
for (j in 1:n) {
    jp.tdf[j] <- JP.df(q.4[j], jp.mle$mu, jp.mle$kappa, jp.mle$psi)
}

vm.tqf <- qvonmises(edf(q.4), vm.mle$mu, vm.mle$kappa, from = circular(0), tol = 1e-06)
jp.tqf <- 0
for (j in 1:n) {
    jp.tqf[j] <- JP.qf(edf(q.4)[j], jp.mle$mu, jp.mle$kappa, jp.mle$psi)
}

# P-P plot (will magnify deviations in the centre of the plot)
pdf(file = "PP-plot.pdf")
plot(c(0,1), c(0,1), type = "l", col = "black", cex.axis = 1.6, cex.lab = 1.6, xlab = "Fitted distribution function", ylab = "Empirical distribution function")
points(vm.tdf, edf(q.4), pch = 20, cex = 1.2, col = vM.colour)
points(jp.tdf, edf(q.4), pch = 4, lwd = 2, col = JP.colour)
legend("bottomright", bty = "n", pch = c(20, 4), col = c(vM.colour, JP.colour),
       legend = c("von Mises candidate", "Jones-Pewsey candidate"), cex = 1.6)
dev.off()

# P-P plot (will magnify deviations in the tails of the plot)
pdf(file = "QQ-plot.pdf")
plot(c(0,2*pi), c(0,2*pi), type = "l", xaxt = "none", yaxt = "none", col = "black", cex.axis = 1.6, cex.lab = 1.6, xlab = "Fitted quantile function", ylab = "Empirical quantile function")
points(matrix(vm.tqf), matrix(q.4), pch = 20, cex = 1.2, col = vM.colour)
points(matrix(jp.tqf), matrix(q.4), pch = 4, col = JP.colour)
legend("bottomright", bty = "n", pch = c(20, 4), col = c(vM.colour, JP.colour),
       legend = c("von Mises candidate", "Jones-Pewsey candidate"), cex = 1.6)
axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
axis(2, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
dev.off()

# residual = observed - expected
vm.pp.res <- edf(q.4) - vm.tdf
vm.qq.res <- q.4 - vm.tqf
jp.pp.res <- edf(q.4) - jp.tdf
jp.qq.res <- q.4 - jp.tqf

pp <- as.data.frame(cbind(q = findInterval(edf(q.4), c(0.25, 0.5, 0.75)),
                          vm.pp = vm.pp.res,
                          jp.pp = jp.pp.res))
pp.summ <- ddply(pp, .(q), summarize, vm = mean(vm.pp), jp = mean(jp.pp))

qq <- as.data.frame(cbind(q = findInterval(q.4, c(pi/2, pi, 3 * pi / 2)),
                          vm.qq = vm.qq.res,
                          jp.qq = jp.qq.res))
qq.summ <- ddply(qq, .(q), summarize, vm = mean(vm.qq), jp = mean(jp.qq))

# P-P residual plot
pdf(file = "PP-residuals.pdf")
plot(matrix(vm.tdf[order(q.4)]), vm.pp.res[order(q.4)], ylim = c(min(vm.pp.res), 2 * max(vm.pp.res)), type = "o", col = vM.colour, cex = 1.2, pch = 20, cex.axis = 1.5, cex.lab = 1.5,
     xlab = "Fitted distribution function", ylab = "Residual")
points(matrix(jp.tdf[order(q.4)]), jp.pp.res[order(q.4)], type = "o", pch = 20, cex = 1.2, col = JP.colour)
abline(h = 0, col = "black", lwd = 2)
abline(v = c(0.25, 0.5, 0.75), col = "grey")
text(x = c(0.1, 0.4, 0.65, 0.9), y = rep(0.15, 4), cex = 1.5, col = vM.colour, round(pp.summ$vm, 3))
text(x = c(0.1, 0.4, 0.65, 0.9), y = rep(0.13, 4), cex = 1.5, col = JP.colour, round(pp.summ$jp, 3))
dev.off()

# Q-Q residual plot
pdf(file = "QQ-residuals.pdf")
plot(matrix(vm.tqf[order(q.4)]), vm.qq.res[order(q.4)], type = "o", ylim = c(min(vm.qq.res), 2 * max(vm.qq.res)), pch = 20, col = vM.colour, cex = 1.2, cex.axis = 1.5, cex.lab = 1.5,
     xlab = "Fitted quantile function", ylab = "Residual")
points(matrix(jp.tqf[order(q.4)]), jp.qq.res[order(q.4)], type = "o", pch = 20, col = JP.colour, cex = 1.2)
abline(h = 0, col = "black", lwd = 2)
abline(v = c(pi/2, pi, 3 * pi / 2), col = "grey")
text(x = c(pi/4, 3 * pi / 4, 5 * pi / 4, 7 * pi / 4), y = rep(0.4, 4), cex = 1.5, col = vM.colour, round(qq.summ$vm, 3))
text(x = c(pi/4, 3 * pi / 4, 5 * pi / 4, 7 * pi / 4), y = rep(0.3, 4), cex = 1.5, col = JP.colour, round(qq.summ$jp, 3))
dev.off()


