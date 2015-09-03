# require: spatstat; maptools; rgeos; FNN; jpeg; plyr; shape; raster
#=================================================================================================
# CONVERT JPEG IMAGE TO BLACK-AND-WHITE RASTER DATA
#=================================================================================================
# 3.1.1 import map image from JPEG. Produces a black-and-white raster and a table of feature types.
import.map <- function(jpg.file, plot = T, threshold = 0.2) {
    
    # load JPEG image and convert to matrix, then raster
    jpegimage <- readJPEG(jpg.file)
    data <- matrix(1-jpegimage, nrow = nrow(jpegimage), ncol = ncol(jpegimage))
    r <- raster(data, xmn = 0, ymn = 0, xmx = 1, ymx = nrow(data)/ncol(data))
    
    # convert to binary raster
    r0 <- r >= threshold
    
    # identify clumps of points
    cc <- clump(r0, dir = 8)
    
    # create table to store clump 'type'
    clump.types <- cbind(unique(getValues(cc)), 0)
    clump.types <- clump.types[!is.na(clump.types[,1]),]
    
    if (plot) {
        org.par <- par()
        par(mar = c(0,0,0,0))
        plot(cc, col = "black", asp = T, legend = F)
        par(mar = org.par$mar)
    }
    list(features = cc, feature.types = clump.types)
}

#-------------------------------------------------------------------------------------------------
# save feature raster with associated table of feature types
save.features <- function(obj.to.save, name) {
    writeRaster(obj.to.save$features, 
                filename = paste(name, "-features.grd", sep = ""), overwrite = T)
    write.csv(obj.to.save$feature.types, 
              file = paste(name, "-feature-types.csv", sep = ""), row.names = F)
}

#-------------------------------------------------------------------------------------------------
# load feature raster with associated table of feature types
load.features <- function(name) {
    feature.types <- read.csv(file = paste(name, "-feature-types.csv", sep = ""))
    features <- raster(paste(name, "-features.grd", sep = ""))
    
    list(features = features, feature.types = feature.types)
}


#=================================================================================================
# IDENTIFY SCALE AND N-S MARKER
#=================================================================================================
# 3.1.2: identify scale marker and rescale map accordingly
get.scale <- function(site, l = 17) {
    
    ct <- site$feature.types
    cc <- site$features
    cs <- cbind(ct[,1],0)
    
    # recreate binary data
    r0 <- reclassify(cc, rcl = rbind(cbind(ct[,1], 1), c(NA, 0)))
    
    org.par <- par(); par(mar = c(0,0,0,0))
    
    if (length(ct[ct[,2] == 2, 1]) == 1) {
        # a feature has already been classified as the scale
        candidate <- ct[ct[,2] == 2, 1]
        cs[candidate, 2] <- 2
        
        # show it on the map
        plot(reclassify(cc, rcl = cs), col = c("Black", "Red"), asp = T, legend = F)
        par(mar = org.par$mar)
    } else {
        # filter to look for long horizontal lines of black cells
        scale.filter <- matrix(c(0,1,0), nrow = 3, ncol = l)
        
        # only picks up the highest-scoring regions
        cf <- table(cc[focal(r0, w = scale.filter, pad = TRUE, padValue = 0) == l])
        
        # find the largest of the highest-scoring regions, set to 2
        candidate <- as.numeric(names(cf)[which.max(cf)])
        reset <- ct[as.numeric(names(cf)[which.max(cf)]),2]
        ct[candidate,2] <- 2 
        cs[candidate,2] <- 2
        
        # show it on the map
        plot(reclassify(cc, rcl = cs), col = c("Black", "Red"), asp = T, legend = F)
        par(mar = org.par$mar)
        
        # get confirmation: is this the right scale item?
        conf <- readline(prompt = "Is the highlighted object the map scale? (y/n) > ")
        # if not, user selects scale marker manually by clicking on map
        if (tolower(substr(conf, 1, 1)) != "y") {
            ct[candidate,2] <- reset
            cs[candidate,2] <- 0
            print("Please click on the scale marker...", quote = F)    
            scale.clump <- click(cc, n = 1, xy = T, id = F, show = F)
            
            # If not on a point, get nearest point 
            if (is.na(scale.clump$value)) {
                coords <- cbind(getValues(cc), xyFromCell(cc, cell = 1:ncell(cc)))
                coords <- coords[!(is.na(coords[,1])), ]
                nn <- knnx.index(coords[,2:3], data.frame(scale.clump$x, scale.clump$y), k = 1)
                candidate <- coords[as.numeric(nn),1]
            } else {
                candidate <- scale.clump$value
            }  
            ct[candidate,2] <- 2
            cs[candidate,2] <- 2
            plot(reclassify(cc, rcl = cs), col = c("Black", "Red"), asp = T, legend = F)
        }
    }           
    # get real length of scale
    m.scale <- "x"
    while (is.na(suppressWarnings(as.numeric(m.scale)))) {
        m.scale <- readline(prompt = "What is the length, in metres, represented by the map scale? > ")
    }
    m.scale <- as.numeric(m.scale)
    
    # get length of scale marker feature
    sc <- xyFromCell(cc, 1:ncell(cc))
    sc <- data.frame(sc[(!is.na(getValues(cc))) & (getValues(cc) == candidate),])
    scale.length <- max(sc$x) - min(sc$x)
    sc.check <- min(sc$x) + scale.length
    
    rescale <- readline(prompt = "Do you want to rescale the map coordinates now? > ")
    if (tolower(substr(rescale, 1, 1)) != "y") {
        list(true.length = m.scale, scale.length = scale.length, clump.id = candidate)
    } else {
        scale <- scale.length / m.scale
        cc.xy <- rasterFromXYZ(cbind(xyFromCell(cc, cell = 1:ncell(cc))/ scale,
                                     z = getValues(cc)))
        print("New object 'rescaled' created containing rescaled raster.", quote = F)
        rescaled <<- list(features = cc.xy, feature.types = ct)
    }
}

#-------------------------------------------------------------------------------------------------
# 3.1.2: identify N-S marker and return its angle in radians
get.NS.axis <- function(site, show.result = F) {
    
    cc <- site$features
    ct <- site$feature.types
    coords <- xyFromCell(cc, cell = 1:ncell(cc))
    coords <- cbind(getValues(cc), coords)
    coords <- coords[!(is.na(coords[,1])), ]
    
    if (length(ct[ct[,2] == 3, 1]) == 1) {
        # already identified - just need to extract angle
        target <- ct[ct[,2] == 3, 1]
    } else {
        # manually select object to use as marker
        show.features(site)
        print("Please click on the N-S marker...", quote = F)    
        scale.clump <- click(site$features, n = 1, xy = T, id = F, show = F)
        
        # If not on a point, get nearest point 
        if (is.na(scale.clump$value)) {
            nn <- knnx.index(coords[,2:3], data.frame(scale.clump$x, scale.clump$y), k = 1)
            target <- coords[as.numeric(nn),1]
        } else {
            target <- scale.clump$value
        }
    }
    
    # extract points in selected shape
    coords <- coords[coords[,1] == target, 2:3]
    
    # idenfity longest axis
    L <- coords[which.min(coords[,1]),]; R <- coords[which.max(coords[,1]),]
    B <- coords[which.min(coords[,2]),]; U <- coords[which.max(coords[,2]),]
    
    h <- sqrt(sum((L-R)^2)); v <- sqrt(sum((U-B)^2))
    if (v > h) {
        axis <- atan2((U-B)[2], (U-B)[1])
        x.0 <- B[1]; y.0 <- B[2]; l <- v * 0.75
    } else {
        axis <- atan2((R-L)[2], (R-L)[1])
        x.0 <- L[1]; y.0 <- L[2]; l <- h * 0.75
    }
    # set feature type to 3
    ct[target, 2] <- 3
    
    angle.found <- list(features = site$features, feature.types = ct)
    if (show.result) {
        show.features(angle.found)
        r <- mean(res(site$features)) * 100
        Arrows(x0 = x.0, y0 = y.0, 
               x1 = x.0 + l * cos(axis), y1 = y.0 + l * sin(axis),
               code = 2, lwd = 2, col = "blue")
    }
    # create new object containing updated feature set
    print("New object 'NS.marked' created with axis marked.", quote = F)
    NS.marked <<- angle.found
    axis
}

#=================================================================================================
# IDENTIFY NON-POST-HOLE FEATURES
#=================================================================================================
# 3.1.3: exclude sparse shapes from potential post-hole set
exclude.sparse.shapes <- function(site, density = 0.55, lower = 3, upper, plot = T) {
    
    # find annotations/smaller linear features: long, thin shapes
    # get sizes of all clumps identified
    cc <- site$features; ct <- site$feature.types
    fr <- freq(cc)
    fr <- cbind(fr[!is.na(fr[,1]),], sq.density = 0)
    
    # if no limit set for large features, set limit just above size of largest feature 
    if (missing(upper)) {upper <- max(fr[,2] + 1)}
    
    # get table of x and y coordinates by feature
    xy <- cbind(xyFromCell(cc, cell = 1:ncell(cc)), id = getValues(cc))
    xy <- as.data.frame(xy[!is.na(xy[,3]),])
    
    # get height and width of each feature, in pixels
    for (i in 1:nrow(fr)) {
        coords <- matrix(xy[xy[,3] == i, ], ncol = 3)
        max.l <- max(length(unique(coords[,1])),length(unique(coords[,2])))
        
        # shape density: if shape was bounded by a square, how much is filled?
        fr[i,3] <- round(fr[i,2] / max.l^2,3)
    }        
    # very small 'noise' features and sparse features are classed as annotations
    ct[(ct[,2] <= 1) & (fr[,2] < lower), 2] <- 4
    ct[(ct[,2] <= 1) & (fr[,3] < density), 2] <- 4
    
    # assign large features separately - won't be tested for neighbours later
    ct[(ct[,2] %in% c(0, 4)) & (fr[,2] >= upper), 2] <- 5
    
    # create a new object containing the updated feature types
    sparse.shapes.classified <<- list(features = cc, feature.types = ct)
    if (plot) {show.features(sparse.shapes.classified)
               title (main = "Features reclassified according to size & compactness")}
}

#-------------------------------------------------------------------------------------------------
# 3.1.4: exclude complex shapes from potential post-hole set
feature.closing <- function(site, progress.bar = F, plot.progress = F, rm = 1) {
    fts <- site$features
    
    # select only potential post-hole candidates for conversion
    cand <- cbind(site$feature.types[,1], site$feature.types[,1])
    cand[site$feature.types[,2] > 0, 2] <- NA
    
    # radius of B: 1 cell width, modified by rm if necessary
    r <- signif(mean(res(fts)),2) * rm
    
    # convert feature points to polygons
    site.poly <- rasterToPolygons(reclassify(fts, cand), dissolve = T)
    if (plot.progress) {plot(site.poly)}
    
    dilated.areas <- matrix(nrow = length(site.poly), ncol = 3)
    closed <- list()
    if (progress.bar) {pb <- txtProgressBar(min = 0, max = length(site.poly), style = 3)}
    
    # support functions (because coercion function has known bug)
    # taken from https://stat.ethz.ch/pipermail/r-sig-geo/2009-May/005781.html
    owin2Polygons <- function(x, id="1") {
        stopifnot(is.owin(x))
        x <- as.polygonal(x)
        closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
        pieces <- lapply(x$bdry,
                         function(p) {
                             Polygon(coords=closering(cbind(p$x,p$y)),
                                     hole=is.hole.xypolygon(p))  })
        z <- Polygons(pieces, id)
        return(z)
    }
    
    owin2SP <- function(x) {
        stopifnot(is.owin(x))
        y <- owin2Polygons(x)
        z <- SpatialPolygons(list(y))
        return(z)
    }
    
    # support function to remove holes from polygons
    # taken from https://stat.ethz.ch/pipermail/r-sig-geo/2014-January/020139.html
    remove.polygon.holes <- function(SPDF) {
        # SPDF can be a SpatialPolygonsDataFrame or a SpatialPolygons object
        polys <- slot(SPDF, "polygons")
        holes <- lapply(polys, function(x) sapply(slot(x, "Polygons"), slot,
                                                  "hole"))
        res <- lapply(1:length(polys), function(i) slot(polys[[i]],
                                                        "Polygons")[!holes[[i]]])
        IDs <- row.names(SPDF)
        SpatialPolygons(lapply(1:length(res), function(i)
            Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(SPDF)))
    }
    
    for (i in 1:length(site.poly)) {
        f <- remove.polygon.holes(site.poly[i,])
        closed[[i]] <- closing(as(f, "owin"), r, polygonal = T)
        
        dilated.areas[i,1] <- unique(getValues(fts)[cellFromPolygon(fts, site.poly[i,])[[1]]])
        dilated.areas[i,2] <- length(cellFromPolygon(fts, site.poly[i,])[[1]])
        dilated.areas[i,3] <- length(cellFromPolygon(fts, owin2SP(closed[[i]]))[[1]])
        
        if (plot.progress) {
            # colour plot according to difference in area after closing
            if (abs(dilated.areas[i,2] - dilated.areas[i,3]) > 1) {cl = "red"} else {cl = "blue"}
            plot(closed[[i]], col = adjustcolor(cl, alpha.f = 0.5), add = T)
        }
        if (progress.bar) {setTxtProgressBar(pb, i)}
    }
    
    fts <- site$feature.types
    fts[fts[,1] %in% dilated.areas[dilated.areas[,3] - dilated.areas[,2] > 1,1],2] <- 4
    features.closed <<- list(features = site$features, feature.types = fts)
    
    list(ft.areas = dilated.areas, closing = closed)
}

#-------------------------------------------------------------------------------------------------
# support function to get dimensions of all features in feature raster
feature.dims <- function(site) {
    
    xy <- data.frame(cbind(id = getValues(site$features), 
                           xyFromCell(site$features, 1:ncell(site$features))))
    xy <- xy[!is.na(xy[,1]),]
    
    dims <- ddply(xy, .(id), summarise, xm = mean(x), ym = mean(y),     # mid-points
                  x.max = max(x), x.min = min(x),                       # left and right edge
                  y.max = max(y), y.min = min(y),                       # top and bottom edge
                  freq = length(x))
    
    # add width & height of shape, in cells
    dims <- cbind(dims, width = (dims$x.max - dims$x.min) / res(site$features)[1] + 1,
                  height = (dims$y.max - dims$y.min) / res(site$features)[2] + 1)
    
    d.max <- apply(cbind(dims$width, dims$height), 1, max)
    d.min <- apply(cbind(dims$width, dims$height), 1, min)
    dims$density <- dims$freq/(d.max^2)         # shape density
    dims$ratio <- d.min/d.max                   # ratio of shortest to longest edge
    dims
}

#-------------------------------------------------------------------------------------------------
# 3.1.5: fill in broken site boundary
fill.broken.boundary <- function(site, plot.progress = F, s = 0.2) {
    z <- feature.dims(site)
    if (plot.progress) {plot(reclassify(site$features, cbind(z$id, z$density < s)), legend = F)}
    # look at very sparse features in turn
    sp <- z$id[z$density < s]
    cc <- site$features
    nbs <- c()
    n <- c()
    for (i in 1:length(sp)) {
        xy <- xyFromCell(cc, Which(cc == sp[i], cells = TRUE))
        ext <- extent(xy)
        if (plot.progress) {points(xy, pch = 20, cex = 0.2, col = "black")}
        
        # find ends of lines
        if ((ext@xmax - ext@xmin) > (ext@ymax - ext@ymin)) {
            e1 = xy[which.max(xy[,1]),]; e2 = xy[which.min(xy[,1]),]
        } else {
            e1 = xy[which.max(xy[,2]),]; e2 = xy[which.min(xy[,2]),]
        }
        # get length and angle of feature, use to define transect line
        l <- sqrt(sum((e1 - e2)^2))
        a <- atan2(e1[2] - e2[2], e1[1] - e2[1])
        adj <- c(l * cos(a), l * sin(a))
        ends <- rbind(e1 + adj, e2 - adj)
        transect <- SpatialLines(list(Lines(list(Line(ends)), ID = 1)))
        if (plot.progress) {lines(transect, col = "blue")}
        # record neighbouring features
        tmp <- unique(unlist(extract(cc, transect)))
        nbs <- c(nbs, tmp[!is.na(tmp) & tmp != sp[i]])
        n[i] <- length(tmp[!is.na(tmp) & tmp != sp[i]])
    }
    # boundary contains any features between two or more linear features on list
    # and features with neighbours on the list
    b <- unique(c(as.numeric(rownames(table(nbs))[table(nbs) > 1]), sp[n > 0]))
    
    # remove any unusually large features (outliers)
    l <- (1.5 * IQR(z$freq[b])) + quantile(z$freq[b], 0.75)
    boundary <- z$id[z$id %in% b & z$freq <= l]
    site$feature.types[site$feature.types[,1] %in% boundary, 2] <- 4
    boundary.filled <<- site
}


#=================================================================================================
# CONVERT POST-HOLES TO ANGLES
#=================================================================================================
# 3.2.1: identify post-hole centres
get.postholes <- function(site, plot = T) {
    cc <- site$features
    ct <- site$feature.types
    
    # get mid-points of all clumps
    xy <- data.frame(cbind(id = getValues(cc), xyFromCell(cc, 1:ncell(cc))))
    mids <- ddply(xy[!is.na(xy$id),], .(id), summarise, xm = mean(x), ym = mean(y))
    
    if (length(ct[ct[,2]==1, 2]) == 0) {
        # if no post-holes positively identified, set all unclassified features as post-holes
        ct[ct[,2] == 0, 2] <- 1
        final.classification <<- list(features = cc, feature.types = ct)
    }
    # create a new table of midpoints, from post-holes only
    centres <<- mids[ct[,2] == 1,2:3]
    if (plot) {
        x.lim <- c(floor(min(xy[,2])/10) * 10, ceiling(max(xy[,2])/10) * 10)
        y.lim <- c(floor(min(xy[,3])/10) * 10, ceiling(max(xy[,3])/10) * 10)
        
        plot(cc, col = "grey", asp = T, legend = F, xlim = x.lim, ylim = y.lim)
        points(centres[,1], centres[,2], pch = 20, cex = 0.5, asp = T, xlim = x.lim, ylim = y.lim)
    }   
}

#-------------------------------------------------------------------------------------------------
# 3.2.1: remove remote points
filter.by.distance <- function(pts, plot = F) {
    d <- c(knn.dist(pts, k = 1))                # distance to nearest neighbour
    l <- quantile(d, 0.75) + (1.5 * IQR(d))     # threshold for outliers
    
    remote <- d <= l
    if (plot) {
        plot(pts[remote, ], pch = 20, cex = 0.5, asp = T)
        points(pts[!remote, ], col = adjustcolor("red", alpha.f = 0.5), pch = 20, cex = 0.5)
    }
    remote                                      # return Boolean array to use as filter
}

#-------------------------------------------------------------------------------------------------
# 3.2.2: get angle to k nearest neighbours for each point
k.nearest.angles <- function(pts, k) {
    nn <- knn.index(pts, k = k)         # identify nearest neighbours
    angles <- matrix(NA, nrow = nrow(pts), ncol = k)
    for (i in 1:nrow(pts)) {
        for (j in 1:k) {
            diff.x <- pts[nn[i,j], 1] - pts[i, 1]
            diff.y <- pts[nn[i,j], 2] - pts[i, 2]
            angles[i, j] <- atan2(diff.y, diff.x)
        }
    }
    cbind(pts, angles)
}
