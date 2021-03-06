setwd("~/Documents/ArchStats/Dissertation/sections/data-cleaning/img")
#===========================================================================================
# DIAGRAMS TO ILLUSTRATE DILATION, EROSION & CLOSING
#===========================================================================================

library(maptools)
par(mar = c(2,2,0,0))

r <- mean(res(genlis$features))

xy <- xyFromCell(genlis$features, Which(genlis$features == 99, cell = T))
ext <- extent(xy + 2 * c(-r,r,-r,r))
crop1 <- crop(genlis$features, ext)
poly <- rasterToPolygons(crop1, dissolve = T)
dilated <- dilation(as(poly, "owin"), r, polygonal = T)
closed <- closing(as(poly, "owin"), r, polygonal = T)

pdf("cl-complex-1-org.pdf")
plot(crop1, col = "seagreen", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Original feature", cex.main = 3, xaxt = "none", yaxt = "none")
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("red", alpha.f = 0.5))
plot(as.owin(poly), lwd = 2, add = T, hatch = T)
# add B to plot
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, pch = 20)
text(min(xy[,1]) - (7 * r/4), sort(unique(xy[,2]))[3] + (4 * r/4), "B", cex = 3, font = 3)
dev.off()

pdf("cl-complex-2-dilated.pdf")
plot(crop1, col = "seagreen", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Dilation", cex.main = 3, xaxt = "none", yaxt = "none")
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
plot(dilated, hatch = T, add = T, lwd = 2, col = adjustcolor("skyblue", alpha.f = 0.4))
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("red", alpha.f = 0.5))
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, pch = 20)
text(min(xy[,1]) - (7 * r/4), sort(unique(xy[,2]))[3] + (4 * r/4), "B", cex = 3, font = 3)
dev.off()

pdf("cl-complex-3-eroded.pdf")
plot(crop1, col = "seagreen", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Erosion", cex.main = 3, xaxt = "none", yaxt = "none")
plot(dilated, add = T)
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
plot(closed, hatch = T, add = T, lwd = 2, col = adjustcolor("skyblue", alpha.f = 0.4))
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] + r/2, circles = r, inches = F, add = T, bg = adjustcolor("red", alpha.f = 0.5))
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] + r/2, pch = 20)
text(min(xy[,1]) - (7 * r/4), sort(unique(xy[,2]))[3] + (8 * r/4), "B", cex = 3, font = 3)
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
plot(crop1, col = "seagreen", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Original feature", cex.main = 3, xaxt = "none", yaxt = "none")
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("red", alpha.f = 0.5))
# add B to plot
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, pch = 20)
text(min(xy[,1]) - (7 * r/4), sort(unique(xy[,2]))[3] - (8 * r/4), "B", cex = 3, font = 3)
dev.off()

pdf("cl-simple-2-dilated.pdf")
plot(crop1, col = "seagreen", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Dilation", cex.main = 3, xaxt = "none", yaxt = "none")
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
plot(dilated, hatch = T, add = T, lwd = 2, col = adjustcolor("skyblue", alpha.f = 0.4))
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("red", alpha.f = 0.5))
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[3] - r/2, pch = 20)
text(min(xy[,1]) - (7 * r/4), sort(unique(xy[,2]))[3] - (8 * r/4), "B", cex = 3, font = 3)
dev.off()

pdf("cl-simple-3-eroded.pdf")
plot(crop1, col = "seagreen", asp = T, xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), legend = F, main = "Erosion", cex.main = 3, xaxt = "none", yaxt = "none")
plot(dilated, add = T)
abline(v = seq(from = xmin(ext) - r/2, to = xmax(ext) + r/2, by = r), col = "black", lty = 3)
abline(h = seq(from = ymin(ext) - (3*r/2), to = ymax(ext) + (3*r/2), by = r), col = "black", lty = 3)
plot(closed, hatch = T, add = T, lwd = 2, col = adjustcolor("skyblue", alpha.f = 0.4))
symbols(min(xy[,1]) - r/2, sort(unique(xy[,2]))[2] - r/2, circles = r, inches = F, add = T, bg = adjustcolor("red", alpha.f = 0.5))
points(min(xy[,1]) - r/2, sort(unique(xy[,2]))[2] - r/2, pch = 20)
text(min(xy[,1])- (7 * r/4), sort(unique(xy[,2]))[2] - (8 * r/4), "B", cex = 3, font = 3)
dev.off()