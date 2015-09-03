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
points(centres, pch = 20)
#points(centres[!nn.filter,], col  = "purple", lwd = 2, pch = 4)
points(centres[!dist.filter,], col = "red", lwd = 2)
# points(centres[!nn.filter,], pch = 4, col = "red", lwd = 2)
legend("topleft", legend = c("Post-hole", "Excluded by distance filter", "Exluded by angular filter"), bty = "n",
       pch = c(20, 1), col = c("black", "red"), cex = 1.3)
dev.off()
}
#------------------------------------------------------------------------------------
# 87: export results to .csv
{
ests <- cbind(bc = rbind(bc$mu, bc$rho, A1inv(bc$rho), bc$beta2, bc$alpha2, NA),
              vm = rbind(vm$mu, A1(vm$kappa), vm$kappa, NA, NA, NA),
              jp = rbind(jp$mu %% (2*pi), A1(jp$kappa), jp$kappa, NA, NA, jp$psi))
colnames(ests) <- c("est", "lower", "upper", 
                    "est.vm", "lower.vm", "upper.vm",
                    "est.jp", "lower.jp", "upper.jp")
rownames(ests) <- c("mu", "rho", "kappa", "beta2", "alpha2", "psi")
ests <- round(ests, 3)
write.csv(ests, file = "../../csv/Genlis-ests.csv", row.names = T, quote = T)
}
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

bq <- circular(((cuts[1:bins] + cuts[2:(bins + 1)])/2)[findInterval(q, cuts)])
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
{
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
legend("bottomright", bty = "n", pch = c(20, 4), col = c(vM.colour, JP.colour),
       legend = c("von Mises candidate", "Jones-Pewsey candidate"), cex = 1.6)
dev.off()

# Q-Q plot (will magnify deviations in the tails of the plot)
pdf(file = "QQ-plot.pdf")
plot(c(0,2*pi), c(0,2*pi), type = "l", xaxt = "none", yaxt = "none", col = "black", cex.axis = 1.6, cex.lab = 1.6, xlab = "Fitted quantile function", ylab = "Empirical quantile function")
#abline(v = qq.mean, col = "grey")
#abline(h = qq.mean, col = "grey")
points(matrix(vm.tqf), matrix(q.4), pch = 20, cex = 1.2, col = vM.colour)
points(matrix(jp.tqf), matrix(q.4), pch = 20, col = JP.colour)
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
{
col.a <- "seagreen"; col.b <- "lightseagreen"
    
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
curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, lty = 4, col = JP.colour, lwd = 3)
curve(djonespewsey(x, mu = circular(jp.a$mu), kappa = jp.a$kappa, psi = jp.a$psi), n = 3600, add = T, lty = 1, col = col.a, lwd = 3)
axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
legend("topright", bty = "n", cex = 1.3, col = c("black", JP.colour, col.a), lty = c(1,4,1), lwd = 3,
       legend = c("Kernel density estimate", "Global Jones-Pewsey distribution", "Jones-Pewsey for this quadrant"))
dev.off()

pdf(file = "quad-B-hist.pdf")
hist(matrix(b)[quadrant == 1], xaxt = "none", col = point.colour, breaks = 40, cex.axis = 1.5, ylim = c(0,1),
     xlab = "Transformed angle (radians)", border = "darkgrey", cex.lab = 1.5, main = "", freq = F)
lines(kd.b, col = "black", lwd = 3)
curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, lty = 4, col = JP.colour, lwd = 3)
curve(djonespewsey(x, mu = circular(jp.b$mu), kappa = jp.b$kappa, psi = jp.b$psi), n = 3600, add = T, lty = 1, col = col.b, lwd = 3)
axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
legend("topright", bty = "n", cex = 1.3, col = c("black", JP.colour, col.a), lty = c(1,4,1), lwd = 3,
       legend = c("Kernel density estimate", "Global Jones-Pewsey distribution", "Jones-Pewsey for this quadrant"))
dev.off()
}

#===========================================================================================
# E-M CLUSTERING
#===========================================================================================

mcol1 <- "seagreen"; mcol2 = "skyblue"

x <- circular(seq(0, 2 * pi, 0.01))
components <- matrix(nrow = em.u.vm$k, ncol = length(x))
for (i in 1:em.u.vm$k) {
    components[i, ] <- dvonmises(x, circular(em.u.vm$mu[i]), 
                                 em.u.vm$kappa[i]) * em.u.vm$alpha[i]
}
y.max <- max(c(colSums(components), hist(matrix(b), plot = F, breaks = 40)$density)) * 1.1
labl <- paste("Mixture of", em.u.vm$k, "von Mises")

pdf("mixt-uvm-plot.pdf")
hist(matrix(b), freq = F, ylim = c(0, y.max), main = "", xaxt = "none",
     col = "lightgrey", border = "darkgrey", xlab = labl, xlim = c(0, 2 * pi), 
     breaks = 40, cex.lab = 1.5)
lines(matrix(x), components[1,], col = mcol1, lwd = 3, lty = 3)
lines(matrix(x), components[2,], col = mcol1, lwd = 3, lty = 3)
lines(matrix(x), colSums(components), lwd = 3)

curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, lty = 4, col = JP.colour, lwd = 3)

axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
legend("topright", bty = "n", cex = 1.3, col = c("black", mcol1, JP.colour), lty = c(1,3,4), lwd = 3,
       legend = c("Uniform-von Mises mixture model", "Uniform-von Mises components", "Jones-Pewsey model"))
dev.off()


vm.components <- matrix(nrow = em.vm$k, ncol = length(x))
for (i in 1:em.vm$k) {
    vm.components[i, ] <- dvonmises(x, circular(em.vm$mu[i]), 
                                 em.vm$kappa[i]) * em.vm$alpha[i]
}
pdf("mixt-vm-plot.pdf")
hist(matrix(b), freq = F, ylim = c(0, y.max), main = "", xaxt = "none",
     col = "lightgrey", border = "darkgrey", xlab = labl, xlim = c(0, 2 * pi), 
     breaks = 40, cex.lab = 1.5)
lines(matrix(x), vm.components[1,], col = mcol2, lwd = 3, lty = 3)
lines(matrix(x), vm.components[2,], col = mcol2, lwd = 3, lty = 3)
lines(matrix(x), colSums(vm.components), lwd = 3)

curve(djonespewsey(x, mu = circular(jp.mle$mu), kappa = jp.mle$kappa, psi = jp.mle$psi), n = 3600, add = T, lty = 4, col = JP.colour, lwd = 3)

axis(1, at = c(0,0.5,1,1.5,2) * pi, cex.axis = 1.5, 
     labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                expression(paste("3", pi, "/2")), expression(paste("2", pi))))
legend("topright", bty = "n", cex = 1.3, col = c("black", mcol2, JP.colour), lty = c(1,3,4), lwd = 3,
       legend = c("von Mises mixture model", "von Mises components", "Jones-Pewsey model"))
dev.off()


l <- matrix(nrow = length(mu), ncol = length(data))
m <- l
for (i in 1:length(mu)) {
    l[i, ] <- alpha[i] * pvonmises(data, em.vm$mu[i], em.vm$kappa[i], from = circular(0), tol = 1e-06)
    m[i, ] <- alpha[i] * pvonmises(data, em.u.vm$mu[i], em.u.vm$kappa[i], from = circular(0), tol = 1e-06) 
}
uvm.tdf <- colSums(l)
mvm.tdf <- colSums(m)

pdf("mvm-PP.pdf")
plot(c(0,1), c(0,1), type = "l", col = "black", cex.axis = 1.6, cex.lab = 1.6, xlab = "Fitted distribution function", ylab = "Empirical distribution function")
points(jp.tdf, edf(q.4), pch = 20, lwd = 2, col = JP.colour)
points(uvm.tdf, edf(q.4), pch = 20, lwd = 2, col = mcol1)
points(mvm.tdf, edf(q.4), pch = 20, lwd = 2, col = mcol2)
legend("bottomright", bty = "n", pch = c(20, 4), col = c(mcol1, mcol2, JP.colour),
       legend = c("uniform-von Mises mixture", "von Mises mixture", "Jones-Pewsey"), cex = 1.6)
dev.off()

#===========================================================================================

# plot resulting clustering
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

#===========================================================================================
# apply spatial clustering & create 'heatmap' of clusters

plot(pts[em.clusts == 1,], pch = 4, col = "blue", asp = T)
points(pts[em.clusts == 2,], pch = 1, col = "black")

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
