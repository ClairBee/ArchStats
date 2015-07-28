
setwd("~/Documents/ArchStats/Mini-essays/Jpeg-to-points")
library(AS.preprocessing); library(AS.angles); library(AS.circular); library(shape)

catholme <- load.features("Catholme")


classify.sparse.shapes(catholme, density = ((pi/4) + (1/3))/2, upper = 100, plot = T)
catholme <- sparse.shapes.classified; remove(sparse.shapes.classified)
# a number of map features also identified, but only linear ones. Same cutoff works well here - although minimal text to get in the way on this site.


# has very little text to remove: only scale text, outside of the plot.


# a more efficient method is to go straight to removing tall things.
remove.tall.features(catholme, plot = T)
catholme <- tall.features.removed; remove(tall.features.removed)


remove.isolated(catholme, plot = T)
catholme <- density.filtered; remove(density.filtered)


get.postholes(catholme)
catholme <- final.classification; remove(final.classification)
save.features(catholme, "Catholme-final")


#==============================================================================

#k-distances plot
plot(sort(knn.dist(centres, k = 5)), pch = 20); abline(h = 4.65, col = "blue")
db <- dbscan(centres, eps = 4.65, MinPts = 5)

plot(catholme$features, col = "grey", asp = T)
points(centres[,1], centres[,2], col = c(db$cluster), asp = T, pch = 20, cex = 0.5)
c(0:41)[table(db$cluster) > 25]
