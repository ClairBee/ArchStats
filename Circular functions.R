
library(circular)

#============================================================================
# TESTS OF UNIFORMITY & SYMMETRY
#============================================================================

# Graphical assessment of symmetric modality, using the median
median.symmetry.plot <- function(data, enhance.ends = F) {
    
    # convert to radians
    data <- conversion.circular(data, units = "radians")
    
    median <- median.circular(data)
    n <- length(data)
    
    z <- sort((data - median) / 2, decreasing = F)    
    if (enhance.ends) {z <- z/2}
    
    x <- c(); y <- c()
    
    for (i in 1:n) {
        x[i] <- sin(z[i])
        y[i] <- sin((z[i] + z[n + 1 - i])/2)
    }
    
    plot(x, y, pch = 20, main = "Symmetry plot", xlab = "Quantile", ylab = "Symmetrised quantile")
    abline(h = 0, col = "lightseagreen")
}

# Graphical assessment of uniformity (uniform Q-Q plot)
uniformity.plot <- function(data, extend = F, title = T) {
    if (circularp(data)$units == "degrees") {
        data <- conversion.circular(data, units = "radians")
    }
    data <- data %% (2 * pi)
    n <- length(data)
    
    x <- c(1:n)/(n+1)
    y <- sort(data) / (2 * pi)
    
    if(title) {
        m <- "Uniformity plot"
        s <- 1
        if (extend) {
            floor(n/5)
            x.end <- x[1:floor(n/5)] + 1
            y.end <- y[1:floor(n/5)] + 1
            
            x.start <- x[(n - floor(n/5)):n] - 1
            y.start <- y[(n - floor(n/5)):n] - 1
            
            x <- c(x.start, x, x.end)
            y <- c(y.start, y, y.end)
            m <- "Extended uniformity plot"
        }} else {
            m <- ""
            s <- 0
        }    
    
    par(mar = c(4,4,3,1), mai=c(1.1, 1.1, s, 0))
    plot(x, y, pch = 20, asp = T, xlab = "", ylab = "")
    title(m, cex.main = set.cex.main,
          xlab = "Uniform quantiles", ylab = "Sample quantiles", cex.lab = set.cex.lab)
    abline(a = 0, b = 1, col = "lightseagreen")
}

# Run various tests of uniformity (5% significance level hard-coded)
uniformity.tests <- function(data) {
    list(kuiper.stat = eval(kuiper.test(data, alpha = 0.05))$statistic, 
         kuiper = eval(kuiper.test(data, alpha = 0.05))$statistic > 1.747,
         watson.statistic = eval(watson.test(data, alpha = 0.05))$statistic,
         watson = eval(watson.test(data, alpha = 0.05))$statistic > 0.187,
         rao.statistic = eval(rao.spacing.test(data, alpha = 0.05))$statistic,
         rao = eval(rao.spacing.test(data, alpha = 0.05))$statistic > 140.57,
         rayleigh.s = rayleigh.test(data)$statistic,
         rayleigh.p = round(rayleigh.test(data)$p.val, 3),
         rayleigh = rayleigh.test(data)$p.val < 0.05)
}

# ------------------------------------------------- BOOTSTRAP TEST OF UNIFORMITY

# Reflective symmetry: large-sample test
r.symm.test.stat <- function(data) {
    n <- length(data)
    bar <- get.moments(data)
    
    var <- ((1-bar$a4)/2-(2*bar$a2)+(2*bar$a2/bar$r) * (bar$a3 + (bar$a2*(1-bar$a2)/bar$r)))/n
    ts <- abs(bar$b2/sqrt(var))
    pval <- 2*pnorm(ts, mean = 0, sd = 1, lower = F)
    list(test.statistic = ts, p.val = pval)
}

# Reflective symmetry; bootstrap test
r.symm.test.boot <- function(data, B = 9999, alpha = 0.05) {
    
    n <- length(data)
    absz <- r.symm.test.stat(data)$test.statistic
    
    tbar <- mean(data)
    reflection <- 2*tbar-data
    symm.data <- c(data, reflection)
    
    for (b in 2:(B+1)) { 
        bootstrap.sample <- sample(symm.data, size = n, replace = TRUE)
        absz[b] <- r.symm.test.stat(bootstrap.sample)$test.statistic
    }
    list(p.val = length(absz[absz >= absz[1]]) / (B + 1),
         std.error = qnorm(1-(alpha/2)) * sqrt((p.val*(1-p.val))/(B + 1)))
}

#============================================================================
# VON MISES POINT ESTIMATION
#============================================================================

# support function to extract sample moments
get.moments <- function(data) {
    
    t2 <- trigonometric.moment(data, p = 2, center = T)
    t3 <- trigonometric.moment(data, p = 3, center = T)
    t4 <- trigonometric.moment(data, p = 4, center = T)
    
    list(mu = as.numeric(mean(data)),
         r = rho.circular(data),
         b2 = t2$sin,
         b3 = t3$sin,
         b4 = t4$sin,
         a2 = t2$cos,
         a3 = t3$cos,
         a4 = t4$cos)
}

# bias-corrected point estimates (large-sample asymptotic theory)
bc.point.estimates <- function(data, sig = 0.05, symmetric = F) {
    
    n <- length(data)
    bar <- get.moments(data)
    bar$mu <- bar$mu %% (2*pi)
    r2 <- bar$r^2 ; r4 <- r2^2
    
    qval <- qnorm(1-sig/2)
    
    rho.bc <- bar$r - ((1-bar$a2)/(4*n*bar$r))
    r.SE <- sqrt((1-2*r2+bar$a2)/(2*n))    
    rho.upper <- rho.bc + qval*r.SE
    rho.lower <- rho.bc - qval*r.SE
    
    rho <- c(estimate = rho.bc, lower = rho.lower, upper = rho.upper)
    
    if (symmetric) {
        beta2 <- c(estimate = 0, lower = 0, upper = 0)
    } else {    
        beta2.bc <- bar$b2 + ((bar$b3/bar$r)+(bar$b2/r2)-(2*bar$a2*bar$b2/r4))/n
        beta2.SE <- sqrt((((1-bar$a4)/2)-(2*bar$a2)-(bar$b2^2)+(2*bar$a2/bar$r)*(bar$a3+(bar$a2*(1-bar$a2)/bar$r)))/n)
        beta2.upper <- beta2.bc + qval*beta2.SE
        beta2.lower <- beta2.bc - qval*beta2.SE
        
        beta2 <- c(estimate = beta2.bc, lower = beta2.lower, upper = beta2.upper)
    }
    
    div <- 2*n*r2 
    
    mu.bc <- bar$mu + (bar$b2/div)
    mu.SE <- sqrt((1-bar$a2)/div)
    
    mu.upper <- mu.bc + qval * mu.SE
    mu.lower <- mu.bc - qval * mu.SE
    
    mu <- c(estimate = mu.bc, lower = mu.lower, upper = mu.upper)
    
    alpha2.bc <- bar$a2 - (1-(bar$a3/bar$r)-((bar$a2*(1-bar$a2)+bar$b2*bar$b2)/r2))/n
    alpha2.SE <- sqrt((((1+bar$a4)/2)-(bar$a2*bar$a2)+(2*bar$b2/bar$r)*(bar$b3+(bar$b2*(1-bar$a2)/bar$r)))/n)
    alpha2.upper <- alpha2.bc + qval*alpha2.SE
    alpha2.lower <- alpha2.bc - qval*alpha2.SE
    
    alpha2 <- c(estimate = alpha2.bc, lower = alpha2.lower, upper = alpha2.upper)
    
    list(sig = sig, symmetric = symmetric, mu = mu, rho = rho, beta2 = beta2, alpha2 = alpha2)
}

# bias-corrected point estimates (bootstrap method - extra uncertainty)
bc.boot.point.ests <- function(data, symmetric = F, alpha = 0.05, B = 9999) {
    
    n <- length(data)
    
    bc.ests <- function(data, symmetric = F, n) {
        
        bar <- get.moments(data)    
        
        rho.est <- bar$r - ((1-bar$a2)/(4*n*bar$r))
        r2 <- bar$r^2
        r4 <- r2^2
        
        if (symmetric) {
            beta2.est <- 0
        } else {
            beta2.est <- bar$b2 + ((bar$b3/bar$r)+(bar$b2/r2)-(2*bar$a2*bar$b2/r4))/n
        }
        
        div <- 2*n*r2
        mu.est <- (bar$mu + (bar$b2/div)) %% (2*pi)
        
        alpha2.est <- bar$a2 - (1-(bar$a3/bar$r)-((bar$a2*(1-bar$a2)+bar$b2*bar$b2)/r2))/n
        
        list(mu = mu.est, rho = rho.est, beta2 = beta2.est, alpha2 = alpha2.est)   
    }
    
    #    bar <- get.moments(data)  
    ests <- bc.ests(data, symmetric = symmetric, n = n)
    
    mu.est <- ests$mu
    rho.est <- ests$rho
    beta2.est <- ests$beta2
    alpha2.est <- ests$alpha2
    
    if (symmetric) {
        reflection <- 2 * mu.est - data
        sample.data <- c(data, reflection)
    } else {
        sample.data <- data
    }
    
    for (b in 2 : (B+1)) { 
        bootstrap.sample <- sample(sample.data, size = n, replace = TRUE)
        ests <- bc.ests(bootstrap.sample, symmetric, n)
        
        mu.est[b] <- ests$mu
        rho.est[b] <- ests$rho
        beta2.est[b] <- ests$beta2
        alpha2.est[b] <- ests$alpha2
    }
    
    dist <- 0
    
    if (symmetric) {
        dist <- pi - abs(pi - abs(mu.est - mu.est[1]))        
        sdist <- sort(dist)
        
        mu.lower <- mu.est[1] - sdist[(B+1)*(1-alpha)]
        mu.upper <- mu.est[1] + sdist[(B+1)*(1-alpha)]
    } else {
        if (mu.est[1] < pi) {
            ref <- mu.est[1] + pi
            for (b in 1:(B+1)) {
                dist[b] <- -(pi-abs(pi-abs(mu.est[b]-mu.est[1])))
                if (mu.est[b] > mu.est[1]) {
                    if (mu.est[b] < ref) {
                        dist[b] <- -dist[b]
                    }
                }
            }
        } else {    #  if (mu.est[1] >= pi)
            ref <- mu.est[1] - pi
            for (b in 1:(B+1)) { 
                dist[b] <- pi-abs(pi-abs(mu.est[b]-mu.est[1]))
                if (mu.est[b] > ref) {
                    if (mu.est[b] < mu.est[1]) {
                        dist[b] <- -dist[b]
                    }                    
                }                
            }
        }
        sdist <- sort(dist)
        mu.lower <- mu.est[1]+sdist[(B+1)*(alpha/2)]
        mu.upper <- mu.est[1]+sdist[(B+1)*(1-alpha/2)]
        
        sbeta2.est <- sort(beta2.est)
        
        beta2.lower <- sbeta2.est[(B+1)*(alpha/2)]
        beta2.upper <- sbeta2.est[(B+1)*(1-alpha/2)]
        
        beta2 <- c(beta2.est[1], beta2.lower, beta2.upper)
    }    
    
    mu <- c(mu.est[1], mu.lower, mu.upper) 
    
    srho.est <- sort(rho.est)
    
    rho.lower <- srho.est[(B+1)*(alpha/2)] ; rho.upper <- srho.est[(B+1)*(1-alpha/2)]
    
    salpha2.est <- sort(alpha2.est)
    
    alpha2.lower<- salpha2.est[(B+1)*(alpha/2)]
    alpha2.upper <- salpha2.est[(B+1)*(1-alpha/2)]
    
    rho <- c(rho.est[1], rho.lower, rho.upper)
    alpha2 <- c(alpha2.est[1], alpha2.lower, alpha2.upper)
    
    list(mu = mu, rho = rho, beta2 = beta2, alpha2 = alpha2)    
}

#============================================================================
# VON MISES PLOTS & GOODNESS-OF-FIT
#============================================================================

# produce P-P plot against candidate von Mises
vM.PP <- function(data, mu, kappa, title = T)  {
    edf <- ecdf(data)
    tdf <- pvonmises(data, mu, kappa, from=circular(0), tol = 1e-06)
    
    if (title) {
        t = "P-P plot of data vs von Mises"
        s = 1
    } else {
        t = ""
        s = 0
    }
    
    c.par <- par()
    par(mar = c(4,4,3,1), mai=c(1.1, 1.1, s, 0))
    plot.default(tdf, edf(data), pch=20, xlim=c(0,1), ylim=c(0,1), xlab = "", ylab = "")
    title(t, cex.main = set.cex.main, font.main = 1,
          xlab = "von Mises distribution function", ylab = "Empirical distribution function", cex.lab = set.cex.lab)
    lines(c(0,1), c(0,1), lwd=2, col = "lightseagreen")
    par(mar = c.par$mar, mai = c.par$mai)
}

# produce Q-Q plot against candidate von Mises
vM.QQ <- function(data, mu, kappa, title = T)  {
    edf <- ecdf(data)
    tqf <- qvonmises(edf(data), mu, kappa, from=circular(0), tol = 1e-06)
    
    if (title) {
        t = "Q-Q plot of data vs von Mises"
        s = 1
    } else {
        t = ""
        s = 0
    }
    
    c.par <- par()
    par(mar = c(4,4,3,1), mai=c(1.1, 1.1, s, 0))
    
    plot.default(tqf, data, pch=20, xlim=c(0,2*pi), ylim=c(0,2*pi), xlab = "", ylab = "")
    title(t, cex.main = set.cex.main, font.main = 1,
          xlab = "von Mises quantile function", ylab = "Empirical quantile function", cex.lab = set.cex.lab)
    lines(c(0,2*pi), c(0,2*pi), lwd=2, col = "lightseagreen")
    par(mar = c.par$mar, mai = c.par$mai)
}

# bundle of GoF tests for von Mises distribution
# Rao critical value changes with n!
vM.GoF <- function(data, mu, kappa, rao.limit = 140.57) {
    tdf <- pvonmises(data, circular(mu), kappa, from=circular(0), tol = 1e-06)
    cunif <- circular(2*pi*tdf)
    
    list(watson.vM = watson.test(data, alpha = 0.05, dist = "vonmises")$statistic > 0.066,
         kuiper.unif = kuiper.test(cunif, alpha = 0.05)$statistic > 1.747,
         watson.unif = watson.test(cunif, alpha = 0.05)$statistic > 0.187,
         rao.unif = rao.spacing.test(cunif, alpha = 0.05)$statistic > rao.limit,
         rayleigh.unif = rayleigh.test(cunif)$p.val < 0.05,
         rayleigh.unif.pval = round(rayleigh.test(cunif)$p.val, 3))
}

# bundle of bootstrap GoF tests for von Mises distribution
vM.GoF.boot <- function(data, B = 9999) {
    
    n <- length(data)
    vM.mle <- mle.vonmises(data, bias=TRUE)
    mu.0 <- vM.mle$mu
    kappa.0 <- vM.mle$kappa
    
    tdf <- pvonmises(data, mu.0, kappa.0, from=circular(0), tol = 1e-06)
    cunif <- circular(2*pi*tdf)
    
    unif.test.0 <- rep(0,4)
    nxtrm <- rep(1,4)
    
    unif.test.0[1] <- kuiper.test(cunif)$statistic
    unif.test.0[2] <- watson.test(cunif)$statistic
    unif.test.0[3] <- rao.spacing.test(cunif)$statistic
    unif.test.0[4] <- rayleigh.test(cunif)$statistic
    
    for (b in 2:(B+1)) {
        bootstrap.sample <- rvonmises(n, mu.0, kappa.0)
        vM.mle <- mle.vonmises(bootstrap.sample, bias = T)
        mu.1 <- vM.mle$mu
        kappa.1 <- vM.mle$kappa
        
        tdf <- pvonmises(bootstrap.sample, mu.1, kappa.1, from=circular(0), tol = 1e-06)
        cunif <- circular(2*pi*tdf)
        
        nxtrm[1] <- nxtrm[1] + (kuiper.test(cunif)$statistic >= unif.test.0[1])
        nxtrm[2] <- nxtrm[2] + (watson.test(cunif)$statistic >= unif.test.0[2])
        nxtrm[3] <- nxtrm[3] + (rao.spacing.test(cunif)$statistic >= unif.test.0[3])
        nxtrm[4] <- nxtrm[4] + (rayleigh.test(cunif)$statistic >= unif.test.0[4])
    }
    p.val <- nxtrm/(B+1)
    names(p.val) <- c("kuiper", "watson", "rao", "rayleigh")
    p.val
}

#============================================================================
# EVALUATING A JONES-PEWSEY DISTRIBUTION
#============================================================================

# calculate normalising constant of Jones-Pewsey
JP.NCon <- function(kappa, psi){
    if (kappa < 0.001) {ncon <- 1/(2*pi) ; return(ncon)} 
    else {
        eps <- 10*.Machine$double.eps
        if (abs(psi) <= eps) {ncon <- 1/(2*pi*I.0(kappa)) ; return(ncon)}
        
        
        else {
            intgrnd <- function(x){ (cosh(kappa*psi)+sinh(kappa*psi)*cos(x))**(1/psi) }
            
            ncon <- 1/integrate(intgrnd, lower=-pi, upper=pi)$value
            return(ncon) } }
}

# calculate PDF of Jones-Pewsey
JP.pdf <- function(theta, mu, kappa, psi, ncon){
    if (kappa < 0.001) {pdfval <- 1/(2*pi) ; return(pdfval)}
    else {
        eps <- 10*.Machine$double.eps
        if (abs(psi) <= eps) {
            pdfval <- ncon*exp(kappa*cos(theta-mu)) ; return(pdfval) }
        
        
        else { 
            pdfval <- (cosh(kappa*psi)+sinh(kappa*psi)*cos(theta-mu))**(1/psi)
            pdfval <- ncon*pdfval ; return(pdfval) } }
    
}

# calculate density of Jones-Pewsey
JP.df <- function(theta, mu, kappa, psi, ncon) {
    eps <- 10*.Machine$double.eps
    if (theta <= eps) {dfval <- 0 ; return(dfval)}
    
    else 
        if (theta >= 2*pi-eps) {dfval <- 1 ; return(dfval)} else
            if (kappa < 0.001) {dfval <- theta/(2*pi) ; return(dfval)}
    else {
        if (abs(psi) <= eps) {
            vMPDF <- function(x){ ncon*exp(kappa*cos(x-mu)) }
            dfval <- integrate(vMPDF, lower=0, upper=theta)$value
            return(dfval) }
        
        
        else { 
            dfval <- integrate(JP.pdf, mu=mu, kappa=kappa, psi=psi, ncon=ncon, lower=0, upper=theta)$value
            return(dfval) }
    }
}

# calculate quantiles of Jones-Pewsey
JP.qf <- function(u, mu, kappa, psi, ncon) {
    
    eps <- 10*.Machine$double.eps
    if (u <= eps) {theta <- 0 ; return(theta)}
    
    else 
        if (u >= 1-eps) {theta <- 2*pi-eps ; return(theta)} else
            if (kappa < 0.001) {theta <- u*2*pi ; return(theta)}
    else {
        roottol <- .Machine$double.eps**(0.6)
        qzero <- function(x) {
            y <- JP.df(x, mu, kappa, psi, ncon) - u ; return(y) }
        res <- uniroot(qzero, lower=0, upper=2*pi-eps, tol=roottol)
        theta <- res$root ; return(theta) }
}

# simulate values of Jones-Pewsey
JP.sim <- function(n, mu, kappa, psi, ncon) {
    
    fmax <- JP.pdf(mu, mu, kappa, psi, ncon) ; theta <- 0
    for (j in 1:n) {
        stopgo <- 0
        while (stopgo == 0) {
            u1 <- runif(1, 0, 2*pi)
            pdfu1 <- JP.pdf(u1, mu, kappa, psi, ncon)
            u2 <- runif(1, 0, fmax)
            if (u2 <= pdfu1) { theta[j] <- u1 ; stopgo <- 1 }
        } }
    return(theta)
}

#============================================================================
# JONES-PEWSEY POINT ESTIMATION
#============================================================================

# maximum likelihood estimation
mle.jonespewsey <- function(data) {
    n <- length(data)
    s <- sum(sin(data))
    c <- sum(cos(data))
    muvM <- atan2(s,c) %% (2*pi)
    kapvM <- A1inv(sqrt(s*s+c*c)/n)
    
    JPnll <- function(p){
        mu <- p[1] ; kappa <- p[2] ; psi <- p[3] ; parlim <- abs(kappa*psi)
        if (parlim > 10) {
            y <- 9999.0
            return(y)
        } else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(data, mu, kappa, psi, ncon)))
            return (y)
        }
    }
    
    out <- optim(par=c(muvM, kapvM, 0.001), fn=JPnll, gr = NULL, method = "L-BFGS-B", lower = c(muvM-pi, 0, -Inf), upper = c(muvM+pi, Inf, Inf), hessian = TRUE)
  
    list(maxll = -out$value,
         alpha = alpha,
         mu = out$par[1] %% (2*pi),
         kappa = out$par[2],
         psi = out$par[3],
         HessMat = out$hessian)
}

# asymptotic normal-theory confidence intervals
JP.ci.nt <- function(jp.ests, alpha = 0.05) {
    quant <- qnorm(1-alpha/2)
    infmat <- solve(jp.ests$HessMat)
    standerr <- sqrt(diag(infmat))
    
    list(alpha = alpha,
         mu = c(est = jp.ests$mu, lower = jp.ests$mu-(quant*standerr[1]), upper = jp.ests$mu+(quant*standerr[1])),
         kappa = c(est = jp.ests$kappa, lower = jp.ests$kappa-(quant*standerr[2]), upper = jp.ests$kappa+(quant*standerr[2])),
         psi = c(est = jp.ests$psi, lower = jp.ests$psi-(quant*standerr[3]), upper = jp.ests$psi+(quant*standerr[3]))    )
}

# another method, using profile log-likelihood and asymptotic Chi-squared theory,
# is also given by Pewsey, but not used here.

# parametric bootstrap confidence intervals
JP.ci.boot <- function(data, alpha = 0.05, B = 9999) {
    n <- length(data)
    
    ests <- mle.jonespewsey(data)
    mu.est <- ests$mu
    kap.est <- ests$kappa
    psi.est <- ests$psi
    
    ncon <- JP.NCon(kap.est, psi.est)
    
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    for (b in 2:(B+1)) {
        jpdat <- JP.sim(n, mu.est[1], kap.est[1], psi.est[1], ncon)
        ests <- mle.jonespewsey(jpdat) 
        mu.est[b] <- ests$mu
        kap.est[b] <- ests$kappa
        psi.est[b] <- ests$psi
        setTxtProgressBar(pb, b)  
    }
    
    dist <- pi-abs(pi-abs(mu.est-mu.est[1]))
    sdist <- sort(dist)
    
    mu.lower <- mu.est[1]-sdist[(B+1)*(1-alpha)]
    mu.upper <- mu.est[1]+sdist[(B+1)*(1-alpha)]
    
    skap.est <- sort(kap.est) 
    kap.lower <- skap.est[(B+1)*alpha/2]
    kap.upper <- skap.est[(B+1)*(1-alpha/2)]
    
    spsi.est <- sort(psi.est)
    psi.lower <- spsi.est[(B+1)*alpha/2]
    psi.upper <- spsi.est[(B+1)*(1-alpha/2)]
    
    close(pb)
    list(alpha = alpha,
         mu = c(est = mu.est[1], lower = mu.lower, upper = mu.upper),
         kappa = c(est = kap.est[1], lower = kap.lower, upper = kap.upper),
         psi = c(est = psi.est[1], lower = psi.lower, upper = psi.upper),
         HessMat = HessMat)
}

#============================================================================
# JONES-PEWSEY PLOTS & GOODNESS-OF-FIT
#============================================================================

# P-P plot for goodness-of-fit assessment
JP.PP <- function(data, mu, kappa, psi, title = T) {
    n <- length(data)
    ncon <- JP.NCon(kappa,psi) 
    edf <- ecdf(data) 
    tdf <- 0 
    for (j in 1:n) {tdf[j] <- JP.df(data[j], mu, kappa, psi, ncon)}
    
    if (title) {
        t = "P-P plot of data vs Jones-Pewsey"
        s = 1
    } else {
        t = ""
        s = 0
    }
    
    c.par <- par()
    par(mar = c(4,4,3,1), mai=c(1.1, 1.1, s, 0))
    plot.default(tdf, edf(data), pch=20, xlim=c(0,1), ylim=c(0,1), xlab = "", ylab = "")
    title(t, cex.main = set.cex.main, font.main = 1,
          xlab = "Jones-Pewsey distribution function", ylab = "Empirical distribution function", cex.lab = set.cex.lab)
    lines(c(0,1), c(0,1), lwd=2, col = "lightseagreen")
    par(mar = c.par$mar, mai = c.par$mai)
}

# P-P plot for goodness-of-fit assessment
JP.QQ <- function(data, mu, kappa, psi, title = T) {
    n <- length(data)
    ncon <- JP.NCon(kappa, psi) 
    edf <- ecdf(data) 
    tqf <- 0 
    for (j in 1:n) {tqf[j] <- JP.qf(edf(data)[j], mu, kappa, psi, ncon)}
    
    if (title) {
        t = "Q-Q plot of data vs Jones-Pewsey"
        s = 1
    } else {
        t = ""
        s = 0
    }
    
    c.par <- par()
    par(mar = c(4,4,3,1), mai=c(1.1, 1.1, s, 0))
    plot.default(tqf, data, pch=20, xlim=c(0,2*pi), ylim=c(0,2*pi), xlab = "", ylab = "") 
    title(t, cex.main = set.cex.main, font.main = 1,
          xlab = "Jones-Pewsey quantile function", ylab = "Empirical quantile function", cex.lab = set.cex.lab)
    lines(c(0,2*pi), c(0,2*pi), lwd=2, col = "lightseagreen")
    par(mar = c.par$mar, mai = c.par$mai)
}

# LR test against given psi (ie. against von Mises or cardioid) - NT
JP.psi.LR.test <- function(data, psi.0 = 0, alpha = 0.05) {
    
    null.model = "null"
    if (psi.0 == 0) {null.model <- "von Mises"}
    if (psi.0 == 1) {null.model <- "cardioid"}
    if (psi.0 ==-1) {null.model <- "wrapped Cauchy"}
    
    n <- length(data)
    s <- sum(sin(data))
    c <- sum(cos(data))
    mu.vM <- atan2(s,c) %% (2*pi)
    kappa.vM <- A1inv(sqrt(s*s+c*c)/n)
    
    JPnll.psi <- function(p){
        mu <- p[1]
        kappa <- p[2]
        psi <- p[3]
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)}
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(data, mu, kappa, psi, ncon)))
            return (y)
        }
    }
    
    JPnll.psi.0 <- function(p){
        mu <- p[1]
        kappa <- p[2]
        psi <- psi.0 
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)} 
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(data, mu, kappa, psi, ncon))) 
            return(y) 
        }
    }
    
    out <- optim(par=c(mu.vM, kappa.vM, 0), fn=JPnll.psi, gr = NULL, method = "L-BFGS-B", lower = c(mu.vM-pi, 0, -Inf), upper = c(mu.vM+pi, Inf, Inf))
    maxll1 <- -out$value
    muhat1 <- out$par[1]
    kaphat1 <- out$par[2]
    psihat1 <- out$par[3]
    
    out <- optim(par=c(muhat1, kaphat1), fn=JPnll.psi.0, gr = NULL, method = "L-BFGS-B", lower = c(muhat1-pi, 0), upper = c(muhat1+pi, Inf))
    maxll0 <- -out$value 
    muhat0 <- out$par[1] 
    kaphat0 <- out$par[2]
    
    D <- round(-2*(maxll0-maxll1),3)
    
    pval <- pchisq(D, df= 1, lower.tail=F)
    if (pval < alpha) {outcome = "gives"} else {outcome = "does not give"}
    pval <- round(pval, 3)
    
    comment = paste("Jones-Pewsey ", outcome, " an improvement on ", null.model,
                    " model at the ", alpha * 100, "% level.", sep = "")
    
    comparison <- round(rbind(max.ll = c(maxll0, maxll1),
                        mu = c(muhat0, muhat1),
                        kappa = c(kaphat0, kaphat1),
                        psi = c(psi.0, psihat1)),3)
    colnames(comparison) <- c(null.model, "Jones-Pewsey")
    
    list(D = D, p.val = pval, comment = comment, comparison = comparison)
}

# LR test against given psi (ie. against von Mises or cardioid) - bootstrap
JP.psi.LR.boot <- function(data, psi.0 = 0, B = 9999, alpha = 0.05) {
    
    x <- data 
    n <- length(x)
    s <- sum(sin(x))
    c <- sum(cos(x))
    mu.vM <- atan2(s,c) %% (2*pi)
    kappa.vM <- A1inv(sqrt(s*s+c*c)/n)
    
    JPnll.psi <- function(p){
        mu <- p[1]
        kappa <- p[2]
        psi <- p[3] 
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)}
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(x, mu, kappa, psi, ncon))) ; return (y) 
        }
    }
    
    JPnll.psi0 <- function(p){
        mu <- p[1] 
        kappa <- p[2] 
        psi <- psi.0
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)} 
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(x, mu, kappa, psi, ncon))) ; return(y) }
    }
    
    out <- optim(par=c(mu.vM, kappa.vM, 0), fn=JPnll.psi, gr = NULL, method = "L-BFGS-B", lower = c(mu.vM-pi, 0, -Inf), upper = c(mu.vM+pi, Inf, Inf))
    maxll1 <- -out$value 
    muhat1 <- out$par[1]
    kaphat1 <- out$par[2]
    
    out <- optim(par=c(muhat1, kaphat1), fn=JPnll.psi0, gr = NULL, method = "L-BFGS-B", lower = c(muhat1-pi, 0), upper = c(muhat1+pi, Inf))
    maxll0 <- -out$value 
    muhat0 <- out$par[1]
    kaphat0 <- out$par[2]
    
    D <- -2*(maxll0-maxll1) 
    
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    nxtrm <- 1
    ncon <- JP.NCon(kaphat0, psi.0)
    for (j in 2:(B+1)) {
        stopgo <- 0
        while (stopgo == 0) {
            
            
            x <- JP.sim(n, muhat0, kaphat0, psi.0, ncon)
            
            out <- optim(par=c(muhat0, kaphat0, psi.0), fn=JPnll.psi, gr = NULL, method = "L-BFGS-B", lower = c(muhat0-pi, 0, -Inf), upper = c(muhat0+pi, Inf, Inf))
            maxll1 <- -out$value 
            
            out <- optim(par=c(muhat0, kaphat0), fn=JPnll.psi0, gr = NULL, method = "L-BFGS-B", lower = c(muhat0-pi, 0), upper = c(muhat0+pi, Inf))
            maxll0 <- -out$value
            
            D[j] <- -2*(maxll0-maxll1) 
            
            if (D[j] >= 0) {
                if (D[j] < 40) { 
                    if (D[j] >= D[1]) {nxtrm <- nxtrm+1}
                    stopgo <- 1 } }
        } 
        setTxtProgressBar(pb, j)
    }
    close(pb)
    
    pval <- nxtrm/(B+1) 
    if (pval < alpha) {outcome = "gives"} else {outcome = "does not give"}
    pval <- round(pval,3)
    
    null.model = "null"
    if (psi.0 == 0) {null.model <- "von Mises"}
    if (psi.0 == 1) {null.model <- "cardioid"}
    if (psi.0 ==-1) {null.model <- "wrapped Cauchy"}
    
    comment = paste("Jones-Pewsey ", outcome, " an improvement on ", null.model,
                    " model at the ", alpha * 100, "% level.", sep = "")
    
    comparison <- round(rbind(max.ll = c(maxll0, maxll1),
                              mu = c(muhat0, muhat1),
                              kappa = c(kaphat0, kaphat1),
                              psi = c(psi.0, psihat1)), 3)
    colnames(comparison) <- c(null.model, "Jones-Pewsey")
    
    list(p.val = pval, comment = comment, comparison = comparison)
} 

# AIC & BIC comparison of Jones-Pewsey vs alternative
JP.psi.info <- function(data, psi.0 = 0) {
    null.model = "null"
    if (psi.0 == 0) {null.model <- "von Mises"}
    if (psi.0 == 1) {null.model <- "cardioid"}
    if (psi.0 ==-1) {null.model <- "wrapped Cauchy"}
    
    x <- data
    n <- length(x)
    s <- sum(sin(x)) 
    c <- sum(cos(x))
    mu.vM <- atan2(s,c) %% (2 * pi)
    kappa.vM <- A1inv(sqrt(s*s+c*c)/n)
    
    JPnll <- function(p){
        mu <- p[1] 
        kappa <- p[2] 
        psi <- p[3] 
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)}
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(x, mu, kappa, psi, ncon))) 
            return (y)
        }
    }
    
    JPnllpsi0 <- function(p){
        mu <- p[1] 
        kappa <- p[2] 
        psi <- psi.0 
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)} 
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(x, mu, kappa, psi, ncon))) 
            return(y)
        }
    }
    
    out <- optim(par=c(mu.vM, kappa.vM, 0), fn=JPnll, gr = NULL, method = "L-BFGS-B", lower = c(mu.vM-pi, 0, -Inf), upper = c(mu.vM+pi, Inf, Inf))
    maxll1 <- -out$value 
    muhat1 <- out$par[1]
    kaphat1 <- out$par[2] 
    psihat1 <- out$par[3]
    nu <- 3
    AIC1 <- (2*nu)-(2*maxll1)
    BIC1 <- (log(n)*nu)-(2*maxll1) 
    
    out <- optim(par=c(muhat1, kaphat1), fn=JPnllpsi0, gr = NULL, method = "L-BFGS-B", lower = c(muhat1-pi, 0), upper = c(muhat1+pi, Inf))
    maxll0 <- -out$value 
    nu <- 2
    AIC0 <- (2*nu)-(2*maxll0)
    BIC0 <- (log(n)*nu)-(2*maxll0) 
    
    if (AIC0 < AIC1) {
        comment.AIC <- paste("AIC favours ", null.model, " model over full Jones-Pewsey distributon.", sep = "")
    } else {
        comment.AIC <- paste("AIC favours full Jones-Pewsey distribution over ", null.model, " model.", sep = "")
    }
    if (BIC0 < BIC1) {
        comment.BIC <- paste("BIC favours ", null.model, " model over full Jones-Pewsey distributon.", sep = "")
    } else {
        comment.BIC <- paste("BIC favours full Jones-Pewsey distribution over ", null.model, " model.", sep = "")
    }
    comparison <- round(rbind(mu = c(mu.vM, muhat1),
                              kappa = c(kappa.vM, kaphat1),
                              psi = c(psi.0, psihat1),
                              AIC = c(AIC0, AIC1),
                              BIC = c(BIC0, BIC1)),3)
    colnames(comparison) <- c(null.model, "Jones-Pewsey")
    
    list(comparison, comment.AIC, comment.BIC)
}

# parametric bootstrap goodness-of-fit test for Jones-Pewsey
JP.GoF.boot <- function(data, B = 9999) {
    n <- length(data)
    JPmleRes <- mle.jonespewsey(data) 
    muhat0 <- JPmleRes$mu
    kaphat0 <- JPmleRes$kappa
    psihat0 <- JPmleRes$psi
    
    ncon.0 <- JP.NCon(kaphat0, psihat0) 
    tdf <- 0 
    for (j in 1:n) {
        tdf[j] <- JP.df(data[j], muhat0, kaphat0, psihat0, ncon.0)
    }
    
    cunif <- circular(2*pi*tdf)
    nxtrm <- rep(1,4)
    
    unif.test.0 <- c(kuiper.test(cunif)$statistic, 
                     watson.test(cunif)$statistic,
                     rao.spacing.test(cunif)$statistic,
                     rayleigh.test(cunif)$statistic)
    
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    for (b in 2:(B+1)) {
        bootstrap.sample <- JP.sim(n, muhat0, kaphat0, psihat0, ncon.0)
        JPmleRes <- mle.jonespewsey(bootstrap.sample) 
        muhat1 <- JPmleRes$mu
        kaphat1 <- JPmleRes$kappa
        psihat1 <- JPmleRes$psi
        ncon1 <- JPNCon(kaphat1, psihat1) 
        tdf <- 0
        for (j in 1:n) {
            tdf[j] <- JP.df(bootstrap.sample[j], muhat1, kaphat1, psihat1, ncon1)
        }
        cunif <- circular(2*pi*tdf)
        
        nxtrm[1] <- nxtrm[1] + (kuiper.test(cunif)$statistic >= unif.test.0[1])
        
        nxtrm[2] <- nxtrm[2] + (watson.test(cunif)$statistic >= unif.test.0[2])
        nxtrm[3] <- nxtrm[3] + (rao.spacing.test(cunif)$statistic  >= unif.test.0[3])
        nxtrm[4] <- nxtrm[4] + (rayleigh.test(cunif)$statistic >= unif.test.0[4])
        setTxtProgressBar(pb, b)
    }
    pval <- nxtrm/(B+1)
    names(pval) <- c("kuiper", "watson", "rao", "rayleigh")
    close(pb)
    return(pval)
}

# bundle of GoF tests for Jones-Pewsey distribution
# obtain range of p-values
JP.GoF.pvals <- function(data, mu, kappa, psi) {
    n <- length(data)
    ncon <- JP.NCon(kappa, psi)
    tdf <- 0
    
    for (j in 1:n) {tdf[j] <- JP.df(data[j], mu, kappa, psi, ncon)}
    cunif <- circular(2*pi*tdf)
    
    list(kuiper.res = kuiper.test(cunif),
         watson.res = watson.test(cunif),
         rao.res = rao.spacing.test(cunif),
         rayleigh.res = rayleigh.test(cunif))
}

# obtain critical values for given significance limit
JP.GoF.critical <- function(data, mu, kappa, psi, sig = 0.05) {
    n <- length(data)
    ncon <- JP.NCon(kappa, psi)
    tdf <- 0
    
    for (j in 1:n) {tdf[j] <- JP.df(data[j], mu, kappa, psi, ncon)}
    cunif <- circular(2*pi*tdf)
    
    list(kuiper.res = kuiper.test(cunif, alpha = sig),
         watson.res = watson.test(cunif, alpha = sig),
         rao.res = rao.spacing.test(cunif, alpha = sig),
         rayleigh.res = rayleigh.test(cunif))
}

# obtain T/F for 
JP.GoF <- function(data, mu, kappa, psi, sig = 0.05, rao.limit = 140.57) {
    n <- length(data)
    ncon <- JP.NCon(kappa, psi)
    tdf <- 0
    
    for (j in 1:n) {tdf[j] <- JP.df(data[j], mu, kappa, psi, ncon)}
    cunif <- circular(2*pi*tdf)
    
    list(kuiper.unif = kuiper.test(cunif, alpha = sig)$statistic > 1.747,
         watson.unif = watson.test(cunif, alpha = sig)$statistic > 0.187,
         rao.unif = rao.spacing.test(cunif, alpha = sig)$statistic > rao.limit,
         rayleigh.unif = rayleigh.test(cunif)$p.val < sig,
         rayleigh.unif.pval = round(rayleigh.test(cunif)$p.val, 3))
}


#============================================================================
# CHECK FUNCTIONS AGAINST RESULTS GIVEN IN BOOK
#============================================================================

# set up graphical parameters
set.cex.main <- 1; set.cex.lab <- 1


# p88 - testing for reflective symmetry
r.symm.test.stat(circular(wind))                # p = 4.7720 x10^-8
r.symm.test.stat(circular(fisherB1*2*pi/24))    # p = 0.2090

r.symm.test.boot(circular(fisherB6$set1*2*pi/360))   # p = 0.5391


# p93 - bias-corrected point estimates & confidence intervals, large-sample
bc.point.estimates(circular(wind))     
        # mu    :  0.29 ( 0.20,  0.38)
        # rho   :  0.66 ( 0.60,  0.71)
        # beta2 : -0.20 (-0.27, -0.13)
        # alpha2:  0.43 ( 0.34,  0.52)
bc.point.estimates(circular(fisherB1*2*pi/24))
        # mu    :  4.49 ( 4.21, 4.77)
        # rho   :  0.31 ( 0.24, 0.39)
        # beta2 :     0 (    0,    0)
        # alpha2: -0.05 (-0.13, 0.04)

# p96 - bias-corrected point estimates & confidence intervals, bootstrap method
bc.boot.point.ests(circular(fisherB6$set1*2*pi/360))
        # mu    :  3.98 ( 3.41, 4.55)
        # rho   :  0.39 ( 0.23, 0.56)
        # beta2 :     0 (    0,    0)
        # alpha2: -0.21 (-0.50, 0.13)


# p104 - von Mises P-P and Q-Q plots
cdat <- circular(fisherB6$set1*2*pi/360)
vM.mle <- mle.vonmises(cdat)
vM.PP(cdat, vM.mle$mu, vM.mle$kappa)
vM.QQ(cdat, vM.mle$mu, vM.mle$kappa)

      
# p105 - goodness-of-fit tests
# without allowing for parameter estimation:
# now changed to allow Rao limit to be set manually
unlist(vM.GoF(cdat, vM.mle$mu, vM.mle$kappa))
# Kuiper:   ok
# Watson:   ok
# Rao:      rejects incorrectly----------------- critical value may depend on n?
# Rayleigh: ok

# bootstrap version, allowing for parameter estimation:
vM.GoF.boot(cdat)
# kuiper  watson    rao   rayleigh
# 0.1131  0.1156  0.0333   0.2709 


# p108 - Jones-Pewsey maximum likelihood estimation & normal-theory CI
jp.ci.nt(mle.jonespewsey(cdat), alpha = 0.05)
# max.ll: -64.75623
# mu    :  3.969692     (3.454011  , 4.485374  )
# kappa :  1.448023     (0.2417347 , 2.6543104 )
# psi   :  1.088886     (0.07486634, 2.10290618)

# Hessian matrix: 18.082362 -2.16955022 -2.88511287
#                 -2.169550  2.94986445 -0.08841012
#                 -2.885113 -0.08841012  4.26651631

# Jones-Pewsey bootstrap CI
JP.ci.boot(lcdat)
# mu   : ( 3.462306, 4.477079)
# kappa: ( 0.781502, 11.64282)
# psi  : (-2.266002, 2.180046)

# p113 - likelihood ratio test of JP vs vM/cardioid
JP.psi.LR.test(cdat, psi.0 = 1)
# D: 0.03244359                 p.value: 0.8570575
# max.ll -64.772452 -64.756230
# mu       3.956381   3.969703
# kappa    1.433985   1.448059
# psi      1.000000   1.088908

JP.psi.LR.boot(cdat, psi.0 = 1, B = 9999)
# my values:
# p: 0.883
# max.ll  -60.540      -60.518
# mu        3.956        3.970
# kappa     1.434        1.448
# psi       1.000        1.089

JP.psi.info(cdat, psi.0 = 1)
#       JP      WC      vM      cardioid
# AIC   135.51  139.91  137.32  133.54
# BIC   140.58  143.29  140.70  136.92

JP.GoF.critical(cdat, mu = 3.96, kappa = 1.43, psi = 1, alpha = 0.05)

JP.GoF(cdat, mu = 3.96, kappa = 1.43, psi = 1, rao.limit = 153.82)

JP.GoF.pvals(cdat, mu = 3.96, kappa = 1.43, psi = 1)

JP.GoF.boot(cdat, B = 999)
# watson   kuiper   rao      rayleigh
# 0.159    0.232    0.113    0.293
# 0.165    0.239    0.116    0.320 


