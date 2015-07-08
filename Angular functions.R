

library(circular); library(FNN)


#============================================================================
# EXTRACTING ANGLES BETWEEN POINTS
#============================================================================

# measure angles within a certain radius
pointwise.angles.atan2 <- function(pts, radius) {
    if(missing(radius)) {
        x.r <- max(p[,1]) - min(p[,1])
        y.r <- max(p[,2]) - min(p[,2])
        radius <- sqrt(x.r^2 + y.r^2)
    }
    r <- radius^2
    
    angles <- matrix(nrow = nrow(pts), ncol = nrow(pts))
    for (i in 1:nrow(pts)) {
        for (j in 1:nrow(pts)) {
            if (j == i) {angles[i,j] <- NA}
            else {
                diff.x <- pts[j, 1] - pts[i, 1]
                diff.y <- pts[j, 2] - pts[i, 2]
                
                # if distance is greater than specified radius, ignore
                if (diff.x^2 + diff.y^2 <= r) 
                {angles[i,j] <- atan2(diff.y, diff.x)}
                else {angles[i,j] <- NA}
            }
        }
    }
    cbind(pts,angles)
}


# measure angles to nearest k points
k.nearest.angles <- function(pts, k) {
    
    nn <- knn.index(p, k = k)
    angles <- matrix(NA, nrow = nrow(pts), ncol = nrow(pts))
    for (i in 1:nrow(pts)) {
        for (j in 1:k) {
            nn[i,j]
            diff.x <- pts[nn[i,j], 1] - pts[i, 1]
            diff.y <- pts[nn[i,j], 2] - pts[i, 2]
            angles[i, nn[i,j]] <- atan2(diff.y, diff.x)
        }
    }
    cbind(pts, angles)
}


#============================================================================
# GRAPHICAL REPRESENTATION OF LINKS BETWEEN POINTS
#============================================================================

# highlight selected point & those for which angle has been measured around it
show.selected <- function(angles, n) {
    plot(angles[, 1:2], pch = 20, cex = 0.5, asp = T)
    points(angles[n, 1:2], col = "red")
    points(angles[!is.na(angles[n,-c(1,2)]),1:2], col = "blue")
}

# add line segment showing direction to selected points
show.directions <- function(angles, l = 0.01) {
    
    plot(angles[,1:2], col = "lightgrey", pch = 20, cex = 0.5, asp = T)
    
    for (i in 1:nrow(angles)) {
        for (j in 3:nrow(angles)) {
            if (!is.na(angles[i,j])) {
                x.end <- pts[i,1] + (cos(angles[i,j]) * l)
                y.end <- pts[i,2] + (sin(angles[i,j]) * l)
                segments(pts[i,1], pts[i,2], x.end, y.end)
            }
        }
    }
}


#============================================================================

