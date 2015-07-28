
setwd("~/Documents/ArchStats/Mini-essays/Jpeg-to-points")
library(AS.preprocessing); library(AS.angles); library(AS.circular); library(shape)

genlis <- load.features("Genlis")

classify.sparse.shapes(genlis, density = ((pi/4) + (1/3))/2, upper = 100, plot = T)
genlis <- sparse.shapes.classified; remove(sparse.shapes.classified)

# function not working for some reason
# set 'similar' to false when correcting
remove.annotations(genlis, plot = T)
genlis <- annotations.removed; remove(annotations.removed)


remove.tall.features(genlis, plot = T)
genlis <- tall.features.removed; remove(tall.features.removed)


extend.annotations(genlis)
genlis <- with.extensions; remove(with.extensions)

# probably shouldn't do this if using DBscan for clustering (& is duplication of effort)... test difference in results
remove.isolated2(genlis, plot = T)
genlis <- density.filtered; remove(density.filtered)

get.postholes(genlis)
genlis <- final.classification; remove(final.classification)
show.features(genlis); points(centres[,1], centres[,2], pch = 20, cex = 0.5, asp = T)
save.features(genlis, "Genlis-final")

#==============================================================================
pdf(file = "genlis-features.pdf", width = 7, height = (7 * max(centres[,2])/max(centres[,1]))+1)
cols <- c("grey", "grey", "red3", "red", "gold", "lightseagreen", 
          "purple")
desc <- c("Unclassified", "Post-hole", "Scale", "N-S axis", 
          "Annotation", "Large feature", "Changed")
ind <- sort(unique(genlis$feature.types[, 2])) + 1
l.ind <- ind[(ind != 3) & (ind != 4)]
plot(reclassify(genlis$features, rcl = genlis$feature.types), 
     col = cols[ind], asp = T, legend = F, frame.plot = F)
legend("topleft", legend = c(desc[l.ind], "Mid-point"), col = c(cols[l.ind], "black"),
       pch = 20, cex = 0.6, bty = "n")
points(centres, pch = 20, cex = 0.5, asp = T)
dev.off()
#==============================================================================

#k-distances plot

skd <- sort(knn.dist(centres, k = 10)) 
plot(skd, pch = 20)
abline(a = skd[1], b = (skd[2000] - skd[1])/2000, col = "lightseagreen")
abline(a = -2600 * (skd[2760] - skd[2650])/110, b = (skd[2760] - skd[2650])/110, col = "lightseagreen")

db <- dbscan(centres, eps = 5.03, MinPts = 10)

plot(centres[,1], centres[,2], col = c(db$cluster+1), asp = T, pch = 20)
table(db$cluster)
