setwd("~/Documents/ArchStats/Dissertation/sections/CS1-Genlis/img")
# Row numbers refer to Genlis.R
point.colour <- "grey"; JP.colour = "red"; vM.colour = "blue"; bins <- 90; BW = 30

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
# feature extraction
{
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
points(centres, pch = 20, cex = 1.3)
#points(centres[!nn.filter,], col  = "purple", lwd = 2, pch = 4)
points(centres[!dist.filter,], pch = 20, col = "red", cex = 1.3)
# points(centres[!nn.filter,], pch = 4, col = "red", lwd = 2)
legend("topleft", legend = c("Post-hole", "Excluded by distance filter"), bty = "n",
       pch = 20, col = c("black", "red"), cex = 1.3)
dev.off()
}
#------------------------------------------------------------------------------------
# 87: export results to .csv
{
ests <- cbind(bc = rbind(bc$mu, bc$rho, A1inv(bc$rho), bc$beta2, bc$alpha2, NA),
              vm = rbind(vm$mu, A1(vm$kappa), vm$kappa, NA, NA, NA),
              jp = rbind(jp$mu %% (2*pi), A1(jp$kappa), jp$kappa, NA, NA, jp$psi))
colnames(ests) <- c("est", "lower", "upper", 
                    "estvm", "lowervm", "uppervm",
                    "estjp", "lowerjp", "upperjp")
rownames(ests) <- c("$\\mu$", "$\\rho$", "$\\kappa$", "$\\bar{\\beta}_2$", "$\\bar{\\alpha}_2$", "$\\psi$")
ests <- round(ests, 3)
write.csv(ests, file = "Param-ests.csv", row.names = T, quote = F, na = "-")
}

#------------------------------------------------------------------------------------
# 73: circular and linear plots of the transformed data, with MLE distributions
# bin data
cuts <- c(0:bins) * 2 * pi/bins
b <- circular(((cuts[1:bins] + cuts[2:(bins + 1)])/2)[findInterval(q.4, cuts)])
bq <- circular(((cuts[1:bins] + cuts[2:(bins + 1)])/2)[findInterval(q, cuts)])
kd <- cbind(density.circular(q.4, bw = BW)$x,
            density.circular(q.4, bw = BW)$y)
kd2 <- cbind(x = c(kd[,1], kd[,1] + (2*pi)),
             y = rep(kd[,2], 2))
{
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
plot(bq, axes = F, shrink = 2, col = point.colour, stack = T, sep = 0.05, ylim = c(-1.1,0.9))
lines(density.circular(q, bw = BW), col = "black", lwd = 3)
axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.2,
              labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                         expression(paste("3", pi, "/2"))))
Arrows(0.7 * cos(mx + c(0, pi/2, pi, 3*pi/2)), 0.7 * sin(mx + c(0, pi/2, pi, 3*pi/2)),
       0.8 * cos(mx + c(0, pi/2, pi, 3*pi/2)), 0.8 * sin(mx + c(0, pi/2, pi, 3*pi/2)), lty = 2, col = "seagreen")
legend("bottom", bty = "n", cex = 1.3, col = c("black", "seagreen"), lty = c(1,1,4), lwd = c(3, 2),
       legend = c("Kernel density estimate", expression(paste("Modal direction + 0, ", pi, "/2, ", pi, ", 3", pi, "/2"))))
dev.off()
}


#===========================================================================================
# PLOTS OF VON MISES AND JONES-PEWSEY FITS
#===========================================================================================
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
{
qq.mean <- mean.circular(q.4) %% (2*pi)
pp.mean <- mean(c(JP.df(mean.circular(q.4) %% (2*pi), jp.mle$mu, jp.mle$kappa, jp.mle$psi, jp.ncon),
                  pvonmises(mean.circular(q.4) %% (2*pi), vm.mle$mu, vm.mle$kappa, from = circular(0), tol = 1e-06)))

# P-P plot (will magnify deviations in the centre of the plot)
pdf(file = "PP-plot.pdf")
plot(c(0,1), c(0,1), type = "l", col = "black", cex.axis = 1.6, cex.lab = 1.6, xlab = "Fitted distribution function", ylab = "Empirical distribution function")
#abline(h = JP.df(mean.circular(q.4) %% (2*pi), jp.mle$mu, jp.mle$kappa, jp.mle$psi, jp.ncon), col = "grey")
#abline(v = JP.df(mean.circular(q.4) %% (2*pi), jp.mle$mu, jp.mle$kappa, jp.mle$psi, jp.ncon), col = "grey")
points(vm.tdf, edf(q.4), pch = 20, cex = 1.2, col = vM.colour)
points(jp.tdf, edf(q.4), pch = 20, lwd = 2, col = JP.colour)
legend("bottomright", bty = "n", pch = 20, col = c(vM.colour, JP.colour),
       legend = c("von Mises candidate", "Jones-Pewsey candidate"), cex = 1.6)
dev.off()

# Q-Q plot (will magnify deviations in the tails of the plot)
pdf(file = "QQ-plot.pdf")
plot(c(0,2*pi), c(0,2*pi), type = "l", xaxt = "none", yaxt = "none", col = "black", cex.axis = 1.6, cex.lab = 1.6, xlab = "Fitted quantile function", ylab = "Empirical quantile function")
#abline(v = qq.mean, col = "grey")
#abline(h = qq.mean, col = "grey")
points(matrix(vm.tqf), matrix(q.4), pch = 20, cex = 1.2, col = vM.colour)
points(matrix(jp.tqf), matrix(q.4), pch = 20, col = JP.colour)
legend("bottomright", bty = "n", pch = 20, col = c(vM.colour, JP.colour),
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
plot(matrix(vm.tdf[order(q.4)]), vm.pp.res[order(q.4)], col = vM.colour, cex = 1.2, pch = 20, cex.axis = 1.5, cex.lab = 1.5,
     xlab = "Fitted distribution function", ylab = "Residual")
points(matrix(jp.tdf[order(q.4)]), jp.pp.res[order(q.4)], pch = 20, cex = 1.2, col = JP.colour)
abline(h = 0, col = "black", lwd = 2)
dev.off()

# Q-Q residual plot
pdf(file = "QQ-residuals.pdf")
plot(matrix(vm.tqf[order(q.4)]), vm.qq.res[order(q.4)], pch = 20, col = vM.colour, cex = 1.2, cex.axis = 1.5, cex.lab = 1.5,
     xlab = "Fitted quantile function", ylab = "Residual", xaxt = "none")
#abline(v = mean.circular(q.4) %% 2*pi, col = "grey")
points(matrix(jp.tqf[order(q.4)]), jp.qq.res[order(q.4)], pch = 20, col = JP.colour, cex = 1.2)
axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
abline(h = 0, col = "black", lwd = 2)
dev.off()
}

#===========================================================================================
# QUADRANT PLOTS
#===========================================================================================
col.a <- "seagreen"; col.b <- "purple"
{
pdf(file = "phi-quad-plot.pdf")
plot(bq[quadrant == 0], axes = F, shrink = 2, col = col.a, stack = T, sep = 0.05, ylim = c(-1.1,0.9))
points(bq[quadrant == 1], pch = 20, stack = T, sep = 0.05, shrink = 2, col = col.b)
axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.2,
              labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                         expression(paste("3", pi, "/2"))))
legend("bottom", legend = c("Quadrant A", "Quadrant B"), col = c(col.a, col.b), pch = 20, bty = "n", cex = 1.4)
dev.off()

kd.a <- cbind(density.circular(q.4.a, bw = BW)$x, density.circular(q.4.a, bw = BW)$y)
kd.b <- cbind(density.circular(q.4.b, bw = BW)$x, density.circular(q.4.b, bw = BW)$y)

pdf(file = "quad-A-hist.pdf")
hist(matrix(b)[quadrant == 0], xaxt = "none", col = point.colour, breaks = 40, cex.axis = 1.5, ylim = c(0,1),
     xlab = "Transformed angle (radians)", border = "darkgrey", cex.lab = 1.5, main = "", freq = F)
lines(kd.a, col = "black", lwd = 3)
curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, lty = 4, col = JP.colour, lwd = 4)
curve(djonespewsey(x, mu = circular(jp.a$mu), kappa = jp.a$kappa, psi = jp.a$psi), n = 3600, add = T, lty = 1, col = col.a, lwd = 4)
axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
legend("topright", bty = "n", cex = 1.3, col = c("black", JP.colour, col.a), lty = c(1,4,1), lwd = 4,
       legend = c("Kernel density estimate", "Global Jones-Pewsey distribution", "Jones-Pewsey for this quadrant"))
dev.off()

pdf(file = "quad-B-hist.pdf")
hist(matrix(b)[quadrant == 1], xaxt = "none", col = point.colour, breaks = 40, cex.axis = 1.5, ylim = c(0,1),
     xlab = "Transformed angle (radians)", border = "darkgrey", cex.lab = 1.5, main = "", freq = F)
lines(kd.b, col = "black", lwd = 3)
curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, lty = 4, col = JP.colour, lwd = 4)
curve(djonespewsey(x, mu = circular(jp.b$mu), kappa = jp.b$kappa, psi = jp.b$psi), n = 3600, add = T, lty = 1, col = col.b, lwd = 4)
axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(ip, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
legend("topright", bty = "n", cex = 1.3, col = c("black", JP.colour, col.b), lty = c(1,4,1), lwd = 4,
       legend = c("Kernel density estimate", "Global Jones-Pewsey distribution", "Jones-Pewsey for this quadrant"))
dev.off()
}

#===========================================================================================
# E-M CLUSTERING
#===========================================================================================
mcol1 <- "seagreen"; mcol2 = "skyblue"
# plot mixtures
{
x <- circular(seq(0, 2 * pi, 0.01))
components <- matrix(nrow = em.u.vm$k, ncol = length(x))
for (i in 1:em.u.vm$k) {
    components[i, ] <- dvonmises(x, circular(em.u.vm$mu[i]), 
                                 em.u.vm$kappa[i]) * em.u.vm$alpha[i]
}
y.max <- max(c(colSums(components), hist(matrix(b), plot = F, breaks = 40)$density)) * 1.1

pdf("mixt-uvm-plot.pdf")
hist(matrix(b), freq = F, ylim = c(0, y.max + 0.2), main = "", xaxt = "none",
     col = "lightgrey", border = "darkgrey", xlab = "Transformed angle (radians)", xlim = c(0, 2 * pi), 
     breaks = 40, cex.lab = 1.5)
lines(kd, col = "black", lwd = 3)
lines(matrix(x), components[1,], col = mcol1, lwd = 3, lty = 2)
lines(matrix(x), components[2,], col = mcol1, lwd = 3, lty = 2)
lines(matrix(x), colSums(components), lwd = 3, col = "green4")
curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, lty = 4, col = JP.colour, lwd = 3)
axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
legend("topright", bty = "n", cex = 1.3, col = c("black", "blue", mcol1, JP.colour), lty = c(1,1,2,4), lwd = 3,
       legend = c("Kernel density", "Uniform-von Mises mixture model", "Uniform-von Mises components", "Jones-Pewsey model"))
dev.off()

# plot mixture-vM model (no uniform component)
{
vm.components <- matrix(nrow = em.vm$k, ncol = length(x))
for (i in 1:em.vm$k) {
    vm.components[i, ] <- dvonmises(x, circular(em.vm$mu[i]), 
                                 em.vm$kappa[i]) * em.vm$alpha[i]
}
pdf("mixt-vm-plot.pdf")
hist(matrix(b), freq = F, ylim = c(0, y.max + 0.2), main = "", xaxt = "none",
     col = "lightgrey", border = "darkgrey", xlab = "Transformed angle (radians)", xlim = c(0, 2 * pi), 
     breaks = 40, cex.lab = 1.5)
lines(kd, col = "black", lwd = 3)
lines(matrix(x), vm.components[1,], col = mcol2, lwd = 3, lty = 2)
lines(matrix(x), vm.components[2,], col = mcol2, lwd = 3, lty = 2)
lines(matrix(x), colSums(vm.components), lwd = 3, col = "blue")
curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, lty = 4, col = JP.colour, lwd = 3)
axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
legend("topright", bty = "n", cex = 1.3, col = c("black", "blue", mcol2, JP.colour), lty = c(1,1,2,4), lwd = 3,
       legend = c("Kernel density", "von Mises mixture model", "von Mises components", "Jones-Pewsey model"))
dev.off()
}

data <- q.4
mu <- em.u.vm$mu
l <- matrix(nrow = length(mu), ncol = length(data))
m <- l
for (i in 1:length(mu)) {
    l[i, ] <- em.vm$alpha[i] * pvonmises(data, em.vm$mu[i], em.vm$kappa[i], from = circular(0), tol = 1e-06)
    m[i, ] <- em.u.vm$alpha[i] * pvonmises(data, em.u.vm$mu[i], em.u.vm$kappa[i], from = circular(0), tol = 1e-06) 
}
uvm.tdf <- colSums(l)
mvm.tdf <- colSums(m)

pdf("mvm-PP.pdf")
plot(c(0,1), c(0,1), type = "l", col = "black", cex.axis = 1.6, cex.lab = 1.6, xlab = "Fitted distribution function", ylab = "Empirical distribution function")
points(jp.tdf, edf(q.4), pch = 20, lwd = 2, col = JP.colour)
points(uvm.tdf, edf(q.4), pch = 20, lwd = 2, col = mcol1)
#points(mvm.tdf, edf(q.4), pch = 20, lwd = 2, col = mcol2)
legend("bottomright", bty = "n", pch = 20, col = c(mcol1, JP.colour),
       legend = c("uniform-von Mises mixture", "Jones-Pewsey"), cex = 1.6)
dev.off()
}

# residual plots
{
pdf(file = "PP-mixture-residuals.pdf")
plot(matrix(jp.tdf[order(q.4)]), jp.pp.res[order(q.4)], pch = 20, col = JP.colour, cex = 1.2, cex.axis = 1.5, cex.lab = 1.5,
     xlab = "Fitted distribution function", ylab = "Residual", ylim = c(-0.04, 0.025))
points(matrix(uvm.tdf[order(q.4)]), uvm.pp.res[order(q.4)], pch = 20, col = mcol1, cex = 1.2)
#points(matrix(mvm.tdf[order(q.4)]), mvm.pp.res[order(q.4)], pch = 20, col = mcol2, cex = 1.2)
abline(h = 0, col = "black", lwd = 2)
legend("bottomleft", bty = "n", pch = 20, col = c(mcol1, JP.colour),
       legend = c("uniform-von Mises mixture", "Jones-Pewsey"), cex = 1.6)
dev.off()
}#



#===========================================================================================
# back-transform angles
q4.rep <- rep(q.4/4, 4) + sort(rep(c(0,pi/2, pi, 3*pi/2), length(q.4)))
mvM.rep <- circular(rep((em.u.vm$alpha[2] * dvonmises(q.4, mu = em.u.vm$mu[2], kappa = em.u.vm$kappa[2])) + (em.u.vm$alpha[1] / (2*pi)),4))

# circular plot of the resulting clustering, translated back into phi
{
    pdf(file = "Q-cluster-plot.pdf")
    plot(bq[em.clusts == 1], axes = F, shrink = 2, col = point.colour, stack = T, sep = 0.05, ylim = c(-1.1,0.9))
    points(bq[em.clusts == 2], shrink = 2, col = "black", stack = T, sep = 0.05)
    #lines(density.circular(q, bw = BW), col = "black", lwd = 3)
    #lines.circular(q4.rep, mvM.rep, lwd = 3, col = "blue")
    axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.2,
                  labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                             expression(paste("3", pi, "/2"))))
    legend("bottom", bty = "n", cex = 1.3, col = c("black", "grey"), pch = 20,
           legend = c("von Mises component", "Uniform (noise) component"))
    dev.off()
}

# plot post-holes by cluster
{
    pdf("Genlis-clustered-postholes-1.pdf", height = pdfheight, width = pdfwidth)
    plot(pts, col = c("white", "blue", "green", "red", "orange", "cyan")[db.clust + 1], cex = 1.3, lwd = 3)
    points(pts[em.clusts == 1,], pch = 20, cex = 1.3,  col = "white", xlab = "", ylab = "", cex.axis = 1.3)
    points(pts[em.clusts == 2,], pch = 20, cex = 1.3, col = "black")
    legend("topright", pch = 21, col = c("red", "red", NA, NA), pt.bg = c("black", "white", "black", "darkgrey"), cex = 1.3, pt.cex = 1.3, pt.lwd = 3,
           legend = c("von Mises component in DB cluster", "Uniform component in DB cluster","von Mises component not in DB cluster", "Uniform (noise) component not in DB cluster"), bty = "n")
    points(pts[em.clusts == 1 & db.clust == 0,], cex = 1.3, col = "darkgrey", pch = 20)
    dev.off()
}


#===========================================================================================
# plot the von Mises component
mu <- circular(em.u.vm$mu[2]); kappa <- em.u.vm$kappa[2]; alpha <- em.u.vm$alpha[2]

pdf("fitted-vM-component.pdf")
curve(dvonmises(circular(x), mu, kappa), n = 3600, col = "blue", cex.axis = 1.3, lwd = 3, ylim = c(0,1), xlim = c(0, 2*pi),
      xaxt = "none", xlab = "", ylab = "")
curve(dnorm(x, mu, sqrt(1/kappa)), n = 3600, col = "red", add = T, lwd = 3, lty = 5)
axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
abline(v = qvonmises(c(0.125, 0.875), mu, kappa, from = circular(0)), lwd = 2, lty = 5)   # 75%
abline(v = qvonmises(c(0.95, 0.05), mu, kappa, from = circular(0)), lwd = 2, lty = 4)     # 90%
abline(v = qvonmises(c(0.975, 0.025), mu, kappa, from = circular(0)), lwd = 2, lty = 3)   # 95%
legend("topleft", legend = c("Von Mises", "Normal approximation", "75% of density", "90% of density", "95% of density"),
       lwd = c(3,3,2,2,2), col = c("blue", "red", "black", "black", "black"), lty = c(1,5, 5,4,3), bty = "n", cex = 1.3)

# proportion of measurements that will lie within a particular range....
((2*pvonmises(circular(pi/2), mu = circular(0), kappa)) - 1) * 100      # 99.8% within pi/2 (90 = 22.5)
((2*pvonmises(circular(pi/4), mu = circular(0), kappa)) - 1) * 100      # 91.1% within pi/4 (45 = 11.25)
((2*pvonmises(circular(pi/8), mu = circular(0), kappa)) - 1) * 100      # 61.5% within pi/8 (22.5 = 5.625)



#===========================================================================================

# risk-plot of resulting clustering
{
pts <- centres[dist.filter,]

plot(pts[em.clusts == 2,], pch = 20, col = "black", asp = T)
points(pts[em.clusts == 1,], pch = 20, cex = 0.5, col = "grey", asp = T)

h <- bw.diggle(ppp(x = pts[],1], y = pts[],2], window =  owin(xrange = c(xmin(genlis$features), xmax(genlis$features)),
                                                              yrange = c(ymin(genlis$features), ymax(genlis$features)))))

p.intensity <- function(x, y, obs) {
    ab <- (x - obs[,1])^2 + (y - obs[,2])^2
    (3 / (pi * h^2)) * sum((1 - (ab[ab <= h^2] / h^2))^2)
}

intensity.contour <- function(data) {
    x.plot <- c(floor(xmin(genlis$features)): ceiling(xmax(genlis$features)))
    y.plot <- c(floor(ymin(genlis$features)): ceiling(ymax(genlis$features)))
    z.plot <- matrix(nrow = length(x.plot), ncol = length(y.plot))
    for (i in 1:length(x.plot)) {
        for (j in 1:length(y.plot)) {
            z.plot[i,j] <- p.intensity(x.plot[i], y.plot[j], data)
        }
    }
    contour(x.plot, y.plot, z.plot)
}
intensity.contour(pts[em.clusts == 1,])
intensity.contour(pts[em.clusts == 2,])

risk.contour <- function(points, types, lvls = 5) {
    x.plot <- c(floor(xmin(genlis$features)): ceiling(xmax(genlis$features)))
    y.plot <- c(floor(ymin(genlis$features)): ceiling(ymax(genlis$features)))
    z.plot <- matrix(nrow = length(x.plot), ncol = length(y.plot))
    for (i in 1:length(x.plot)) {
        for (j in 1:length(y.plot)) {
            lam.u = p.intensity(x.plot[i], y.plot[j], points[types == 1,])
            lam.vm = p.intensity(x.plot[i], y.plot[j], points[types == 2,])
            if (lam.vm + lam.u > 0) {
                z.plot[i,j] <- lam.vm / (lam.vm + lam.u)
            } else {
                z.plot[i,j] <- 0
            } } }
    contour(x.plot, y.plot, z.plot, nlevels = lvls)
    list(x = x.plot, y = y.plot, z = z.plot)
}
risk <- risk.contour(pts, em.clusts, lvls = 2)
risk.inv <- risk.contour(pts, 3-em.clusts, lvls = 2)

filled.contour(x = risk$x, y = risk$y, z = risk$z, levels = c(0, .5, 1), col = c("white", "lightgrey"),
               plot.axes = {axis(1); axis(2); 
                            contour(x = risk$x, y = risk$y, z = risk$z, levels = c(0, .5), add = T);
                            points(pts[em.clusts == 1,], pch = 20, col = "red")
                            points(pts[em.clusts == 2,], pch = 20, col = "black")
                            contour(x = risk.inv$x, y = risk.inv$y, z = risk.inv$z, levels = c(0, .5), add = T)})

length(risk$z[risk$z > 0.5])    # 931 / 8468
length(risk.inv$z[risk.inv$z > 0.5])    # 733 / 8468
}
=

#===========================================================================================
# apply spatial clustering & create 'heatmap' of clusters
{
g <- 5
xc <- c(0:ceiling(xmax(genlis$features) / g)) * g
yc <- c(0:ceiling(ymax(genlis$features) / g)) * g

abline(v = xc, col = adjustcolor("seagreen", alpha.f = 0.5))
abline(h = yc, col = adjustcolor("seagreen", alpha.f = 0.5))

qg <- data.frame(count(cbind(x = xc[findInterval(pts[, 1], xc)],
                          y = yc[findInterval(pts[, 2], yc)], z = em.clusts)))
qg.summ <- merge(qg[qg$x.z == 1, c(1,2,4)], qg[qg$x.z == 2, c(1,2,4)], 
                 by = c("x.x", "x.y"), all = T, suffixes = c(".u", ".vm"))
qg.summ[is.na(qg.summ)] <- 0

plot(pts[em.clusts == 1,], pch = 20, col = "grey", asp = T)
points(pts[em.clusts == 2,], pch = 20, col = "black")
points(qg.summ[qg.summ[,4] > 0,] + g/2, pch = 15, col = adjustcolor("red", alpha.f = 0.4), cex = 3)
points(qg.summ[qg.summ[,3] > 0,] + g/2, pch = 15, col = adjustcolor("blue", alpha.f = 0.4), cex = 3)

# proportion of cells containing elements of each cluster
nrow(qg.summ[!is.na(qg.summ[,3]),])/nrow(qg.summ)       # 65% of 5x5 cells contain uniform points
nrow(qg.summ[!is.na(qg.summ[,4]),])/nrow(qg.summ)       # 73% of 5x5 cells contain von Mises points

qg.summ <- cbind(qg.summ, ttl = apply(qg.summ[,3:4],1,sum))
qg.summ <- cbind(qg.summ, prop.u = qg.summ[,3] / qg.summ[,5], prop.vm = qg.summ[,4] / qg.summ[,5])

points(qg.summ[qg.summ[,7] > 0.4, 1:2] + g/2, pch = 15, col = adjustcolor("red", alpha.f = 0.2), cex = 3)
points(qg.summ[qg.summ[,7] > 0.6, 1:2] + g/2, pch = 15, col = adjustcolor("red", alpha.f = 0.2), cex = 3)
points(qg.summ[qg.summ[,7] > 0.8, 1:2] + g/2, pch = 15, col = adjustcolor("red", alpha.f = 0.2), cex = 3)
points(qg.summ[qg.summ[,7] == 1, 1:2] + g/2, pch = 15, col = adjustcolor("red", alpha.f = 0.2), cex = 3)
points(qg.summ[qg.summ[,6] > 0.4, 1:2] + g/2, pch = 15, col = adjustcolor("blue", alpha.f = 0.2), cex = 3)
points(qg.summ[qg.summ[,6] > 0.6, 1:2] + g/2, pch = 15, col = adjustcolor("blue", alpha.f = 0.2), cex = 3)
points(qg.summ[qg.summ[,6] > 0.8, 1:2] + g/2, pch = 15, col = adjustcolor("blue", alpha.f = 0.2), cex = 3)
points(qg.summ[qg.summ[,6] == 1, 1:2] + g/2, pch = 15, col = adjustcolor("blue", alpha.f = 0.2), cex = 3)

# DBscan clustering vs EM clustering: check correlation
library(fpc)
db <- dbscan(pts, MinPts = 4, eps = 4.65)$cluster
plot(sort(knn.dist(pts, k = 4)), pch = 20); abline(h = 4.65)
points(pts[db == 1,], col = "green")
points(pts[db == 2,], col = "cornflowerblue")
points(pts[db == 3,], col = "red")
points(pts[db == 4,], col = "purple")
points(pts[db == 5,], col = "gold")
points(pts[db == 0,], col = "black")
db2 <- db
db2[db2 %in% c(2,4,5,6)] <- 0
plot(pts, pch = 20, asp = T)
points(pts[db2 == 1,], col = "red")
points(pts[db2 == 3,], col = "blue")
cor(db2, em.clusts)


cl <- cbind(genlis$feature.types[,1], NA)
cl[cl[,1] %in% rownames(pts),2] <- em.clusts
r.clusters <- reclassify(genlis$features, cl)

moranL <- MoranLocal(r.clusters)
Geary(r.clusters)

cl[cl[,1] %in% rownames(pts),2] <- db2
Moran(reclassify(genlis$features, cl))
Geary(reclassify(genlis$features, cl))

dir <- cbind(genlis$feature.types[,1], NA)
dir[dir[,1] %in% rownames(pts),2] <- q.4
plot(MoranLocal(reclassify(genlis$features, dir)), legend = F)

lisa <- lisa(pts[,1], pts[,2], em.clusts, neigh = 5)
points(pts[em.clusts == 2,], col = "blue")
}

#===========================================================================================
# proportion of distribution within certain range of pi
pmixt.vonmises(p.to = em.u.vm$mu[2] + (16/180*pi), em.u.vm, p.from = em.u.vm$mu[2])     # 0.1234
pmixt.vonmises(p.to = em.u.vm$mu[2] + (32/180*pi), em.u.vm, p.from = em.u.vm$mu[2])     # 0.2156
pmixt.vonmises(p.to = em.u.vm$mu[2] + (44/180*pi), em.u.vm, p.from = em.u.vm$mu[2])     # 0.2610
pmixt.vonmises(p.to = em.u.vm$mu[2] + (48/180*pi), em.u.vm, p.from = em.u.vm$mu[2])     # 0.2727

#===========================================================================================
# COMPARE DISTRIBUTION OF FIRST AND SECOND NEIGHBOUR ANGLES
{
# same mean, different concentration and therefore distribution. Perhaps unsurprising though.

k.2 <- circular(k.nearest.angles(pts, 2)[,4] %% (2*pi))
q.2 <- (4*k.2) %% (2*pi)

# test for uniformity and symmetry
rayleigh.test(q.2)                  # p = 0         p = 0.0012
kuiper.test(q.2)                    # p < 0.01      p < 0.01
watson.test(q.2)                    # p < 0.01      p < 0.01

r.symm.test.stat(q.2)               # p = 0.374     p = 0.057

bc <- bc.ci.LS(q.2, alpha = 0.05)
# mu     =    3.201099    2.995888    3.406311            3.464893    2.966907    3.962878 
# rho    =   0.3614829   0.2711703   0.4517956          0.16052919  0.06806598  0.25299240
# beta2  = -0.05203799 -0.18135675  0.07728078          0.04734852 -0.10573885  0.20043588
# alpha2 =   0.3003902   0.2135795   0.3872008          0.15547660 -0.02757712  0.33853032 

par(mfrow = c(1,2)); circular.c.plot(q.4); circular.c.plot(q.2)
linear.c.plot(q.4); linear.c.plot(q.2)

q.samples <- list(q.4, q.2)
q.sizes <- c(length(q.4), length(q.2))

watson.common.mean.test(q.samples)                          # p = 0.301
wallraff.concentration.test(q.samples)                      # p = 0.002
mww.common.dist.LS(cs.unif.scores(q.samples), q.sizes)      # p = 0.009
watson.two.test(q.4, q.2)                                   # 0.01 < p < 0.05
watson.two.test.rand(q.4, q.2, NR = 999)                    # p = 0.012
}

# plot of Genlis site with identified grid
pdf("Genlis-grid.pdf", height = pdfheight, width = pdfwidth)
plot(genlis$features, col = "cornflowerblue", cex.axis = 1.3, asp = T, legend = F, frame.plot = F)
points(pts[em.clusts == 2,], pch = 20, cex = 1.3)

hg <- (c(-10:6) * 10)-5
for (i in 1:length(hg)) {
    abline(hg[i], tan(m/4),  col = adjustcolor("gold", alpha.f = 0.5))
}
vg <- (c(1:16) * 10)-5
for (i in 1:length(vg)) {
    abline(vg[i], -tan((m/4)),  col = adjustcolor("gold", alpha.f = 0.5))
}
dev.off()


dev.off()
}