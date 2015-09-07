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
# convert cluster csv into a more useful format
{
csv <- read.csv("DBclust-tests.csv")
csv.unif <- apply(apply(csv[,3:5], 2, "%in%", "Uniform"),1,max)
csv.unif[csv.unif == 0] <- "non-U"; csv.unif[csv.unif == 1] <- "Uniform"
csv.adj <- cbind("size" = csv$size,
                 "unif" = csv.unif,
                 "symm" = levels(csv$symm)[csv$symm],
                 "mu" = round(csv$mu,2),
                 "rho" = round(csv$rho,2),
                 "vm.fit" = apply(cbind(csv$vM.fit.k, csv$vM.fit.w), 1, paste, collapse = ", "),
                 "vm.model" = levels(csv$vM)[csv$vM],
                 "jp.fit" = apply(cbind(csv$JP.fit.k, csv$JP.fit.w), 1, paste, collapse = ", "),
                 "jp.model" = levels(csv$JP)[csv$JP])
csv.adj[csv.adj == "NA, NA" | is.na(csv.adj)] <- "-"
csv.adj <- csv.adj[csv.adj[,9] != "-",c(1,4:9)]
write.table(csv.adj, "DBclust-results.csv", sep = ';', quote = F, row.names = F)
}

#------------------------------------------------------------------------------------
# site plot
cols <- c("blue", "red", "purple", "green2", "lightseagreen")
cat.crop <- crop(catholme$features, c(xmin = 0, xmax = xmax(catholme$features), ymin = 0, ymax = 260))
{
pdf("Catholme-clusters.pdf", height = pdfheight * (260/ymax(catholme$features)), width = pdfwidth)
plot(cat.crop, col = "white", xlim = c(0,250), cex.axis = 1.3, asp = F, legend = F, frame.plot = F)
points(pts, pch = 20, cex = 1.3)
for (i in 1:length(cand)) {
    points(pts[db.clust == cand[i],], col = cols[i], cex = 1.6, pch = 20)
}
legend("bottomleft", ncol = 2, pch = 20, col = cols, legend = paste(l), , pt.cex = 1.6, cex = 1.5)
dev.off()
}
#------------------------------------------------------------------------------------
# plot of cluster mean directions
{
pdf("clust-means.pdf")
plot(circular(0), pch = ".", col = "black", axes = F, shrink = 1.5, xlim = c(-.5,1.5))
axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.2,
              labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                         expression(paste("3", pi, "/2"))))
l <- c()
for (i in 1:length(cand)) {
    m <- JP.mle(q.4[db.clust == cand[i]])$mu
    Arrows(0,0, 0.9*cos(m), 0.9*sin(m), col = cols[i], lwd = 3)
    l[i] <- length(q.4[db.clust == cand[i]])
}
legend("right", ncol = 2, pch = c(rep(20, 5), 1), col = c(cols, "darkgrey"), legend = c(l, "Kernel density estimate"), bty = "n", cex = 1.3) 

dev.off()
}


#------------------------------------------------------------------------------------
# per-cluster circular plot
{
par(mfrow = c(1,6), mar = c(0,0,0,0))
for (i in 1:length(cand)) {
    m <- mean.circular(q.4[db.clust == cand[i]])
    plot(q.4[db.clust == cand[i]], axes = F, shrink = 1.2, cex = 1.3, col = point.colour, stack = T, sep = 0.07)
    axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.2,
                  labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                             expression(paste("3", pi, "/2"))))
    Arrows(0 ,0, 0.9 * cos(m), 0.9 * sin(m), col = cols[i], lwd = 3)
    lines(density.circular(q.4[db.clust == cand[i]], bw = 5), lwd = 3)
}
dev.off()
}

pdf("clust-circular.pdf", width = 6, height = 8)
{
par(mfrow = c(3,2), mar = c(0,0,0,0))
for (i in 1:length(cand)) {
    m <- mean.circular(q.4[db.clust == cand[i]])
    plot(q[db.clust == cand[i]], axes = F, shrink = 1.5, cex = 1.3, col = cols[i], stack = T, sep = 0.07)
    axis.circular(at = c(0,.5,1,1.5) * pi, tcl.text = 0.15, cex = 1.5,
                  labels = c("0", expression(paste(pi, "/2")), expression(paste(pi)),
                             expression(paste("3", pi, "/2"))))
    lines(density.circular(q[db.clust == cand[i]], bw = 5), lwd =  2, col = "darkgrey")
}
plot.new()
legend("center", pch = c(rep(20, 5), NA), lwd = 3, lty = c(rep(NA,5),1),col = c(cols, "darkgrey"), 
       legend = c(l, "Kernel density estimate"), bty = "n", cex = 1.6) 
dev.off()
}
#------------------------------------------------------------------------------------
# plot of angles in each cluster
par(mar = c(2,2,0,0), mfrow = c(2,3))
hist.breaks <- c(20, 20, 20, 20, 20, 5); vm <- c(F,F,T,T,F,T)
{
pdf("clust-models-indiv.pdf", width = 8, height = 6)
par(mar = c(3,3,0,0), mfrow = c(2,3))
for (i in 1:length(cand)) {
    est <- JP.mle(q.4[db.clust == cand[i]])
    hist(matrix(q.4.adj[db.clust == cand[i]]), breaks = hist.breaks[i], col = point.colour, xlim = c(bc$mu[1] + c(-pi,pi)), xaxt = "none",
         ylim = c(0,0.7), cex.axis = 1.5, xlab = "", ylab = "", border = "darkgrey", cex.lab = 1.5, main = "", freq = F)
    axis(1, at = xbreaks, cex.axis = 1.5, labels = xlabl)
    curve(djonespewsey(x, mu = circular(est$mu), kappa = est$kappa, psi = est$psi), n = 3600, add = T, 
          lty = 1, col = cols[i], lwd = 3)
    if (vm[i]) {
        vm.est <- mle.vonmises(q.4[db.clust == cand[i]])
        curve(dvonmises(x, mu = circular(vm.est$mu %% (2*pi)), kappa = vm.est$kappa), n = 3600, add = T, 
              lty = 4, col = cols[i], lwd = 4)
    }
}
    plot.new()
    legend("center", lwd = 3, legend = c("von mises candidate", "Jones-Pewsey candidate"), lty = c(1,4),bty = "n", cex = 1.6) 
dev.off()
}


#=================================================================================================
# FIT ALL POSSIBLE MIXTURE MODELS
#=================================================================================================
# Fit mixture models
write.csv(cbind("cluster" = rownames(clust.AIC), round(clust.AIC,1), signif(sqrt(clust.MSE),2)), "Cluster-results.csv", quote = F, row.names = F)
