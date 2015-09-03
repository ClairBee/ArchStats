setwd("~/Documents/ArchStats/Dissertation/sections/data-cleaning/img")
#===========================================================================================
# DIAGRAMS TO ILLUSTRATE DILATION, EROSION & CLOSING
#===========================================================================================

library(maptools)
select.feature(genlis, replot = T)
par(mar = c(2,2,0,0))

r <- mean(res(genlis$features))

xy <- xyFromCell(genlis$features, Which(genlis$features == 99, cell = T))
ext <- extent(xy + 2 * c(-r,r,-r,r))
crop1 <- crop(genlis$features, ext)
poly <- rasterToPolygons(crop1, dissolve = T)
dilated <- dilation(as(poly, "owin"), r, polygonal = T)
closed <- closing(as(poly, "owin"), r, polygonal = T)

pdf("cl-complex-1-org.pdf")
plot(crop1, col = "grey", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Original feature", cex.main  = 2, xaxt = "none", yaxt = "none")
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("darkgrey", alpha.f = 0.8))
plot(as.owin(poly), lwd = 2, add = T, hatch = T)
# add B to plot
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("darkgrey", alpha.f = 0.8))
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, pch = 20)
dev.off()

pdf("cl-complex-2-dilated.pdf")
plot(crop1, col = "grey", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Dilation", cex.main  = 2, xaxt = "none", yaxt = "none")
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
plot(dilated, hatch = T, add = T, lwd = 2)
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("darkgrey", alpha.f = 0.5))
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, pch = 20)
dev.off()

pdf("cl-complex-3-eroded.pdf")
plot(crop1, col = "grey", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Erosion", cex.main  = 2, xaxt = "none", yaxt = "none")
plot(dilated, add = T)
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
plot(closed, hatch = T, add = T, lwd = 2)
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] + r/2, circles = r, inches = F, add = T, bg = adjustcolor("darkgrey", alpha.f = 0.5))
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] + r/2, pch = 20)
dev.off()

#===========================================================================================
# simpler post-hole feature

xy <- xyFromCell(genlis$features, Which(genlis$features == 360, cell = T))
ext <- extent(xy + 3 * c(-r,r,-r,r))
crop1 <- crop(genlis$features, ext)
poly <- rasterToPolygons(crop1, dissolve = T)
dilated <- dilation(as(poly, "owin"), r, polygonal = T)
closed <- closing(as(poly, "owin"), r, polygonal = T)

pdf("cl-simple-1-org.pdf")
plot(crop1, col = "grey", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Original feature", cex.main  = 2, xaxt = "none", yaxt = "none")
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("darkgrey", alpha.f = 0.8))
plot(as.owin(poly), lwd = 2, add = T, hatch = T)
# add B to plot
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("darkgrey", alpha.f = 0.8))
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, pch = 20)
dev.off()

pdf("cl-simple-2-dilated.pdf")
plot(crop1, col = "grey", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Dilation", cex.main  = 2, xaxt = "none", yaxt = "none")
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
plot(dilated, hatch = T, add = T, lwd = 2)
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("darkgrey", alpha.f = 0.5))
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, pch = 20)
dev.off()

pdf("cl-simple-3-eroded.pdf")
plot(crop1, col = "grey", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Erosion", cex.main  = 2, xaxt = "none", yaxt = "none")
plot(dilated, add = T)
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
plot(closed, hatch = T, add = T, lwd = 2)
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[2] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("darkgrey", alpha.f = 0.5))
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[2] - r/2, pch = 20)
dev.off()