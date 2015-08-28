setwd("~/Documents/ArchStats/Dissertation/img/CS1-Genlis")
# Row numbers refer to Genlis.R
#--------------------------------------------------------------------------------------------
pdfheight <- round(ymax(genlis$features) / 10, 0)
pdfwidth <- round(xmax(genlis$features) / 10, 0)

plot.to.pdf <- function(site, pdf.filename, h, w) {
    site$feature.types[site$feature.types[,2] == 3,2] <- 2
    cols <- c("black", "black", "red", "red", "darkgoldenrod1", "cornflowerblue", "purple")
    desc <- c("Unclassified", "Post-hole", "Scale / N-S axis", "OBSOLETE", "Annotation", "Large feature", "Changed")
    ind <- sort(unique(site$feature.types[, 2])) + 1
    l.ind <- ind[(ind != 3) & (ind != 4)]
    
    pdf(file = pdf.filename, height = h, width = w)
    plot(reclassify(site$features, rcl = site$feature.types), col = cols[ind], cex.axis = 1.3, asp = F, legend = F, frame.plot = F)
    legend("topleft", legend = desc[ind], col = cols[ind], ncol = 2, pch = 20, cex = 1.3, bty = "n", text.font = 3)
    dev.off()
}

# 19: plot sparse shapes classified
plot.to.pdf(sparse.shapes.classified, "Genlis-sparse.pdf", h = pdfheight, w = pdfwidth)

# 23: plot features identified as closed features
plot.to.pdf(features.closed, "Genlis-after-closing.pdf", h = pdfheight, w = pdfwidth)

# 24: plot with broken boundary filled
plot.to.pdf(boundary.filled, "Genlis-boundary-filled.pdf", h = pdfheight, w = pdfwidth)

#28: post-hole centres
pdf("Genlis-1-postholes.pdf", height = pdfheight, width = pdfwidth)
plot(genlis$features, col = "white", cex.axis = 1.3, asp = F, legend = F, frame.plot = F)
points(centres1, pch = 20)
dev.off()

#32: plot annotations removed
plot.to.pdf(annotations.removed, "Genlis-annotations-removed.pdf", h = pdfheight, w = pdfwidth)

#33: tall features removed

#34: annotations extended

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